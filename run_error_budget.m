function budget = run_error_budget(prof_base, imu_ideal, truth, c, cfg_base, N_mc, base_seed)
%RUN_ERROR_BUDGET  Бюджет ошибок: вклад каждого источника в КВО.
%
%   МЕТОД: базовая кампания (все источники включены), затем по одной кампании
%   на каждый источник с ЭТИМ ИСТОЧНИКОМ ОТКЛЮЧЁННЫМ В ГЕНЕРАЦИИ ОШИБОК.
%
%   ПАРНОСТЬ ПО SEED. Все кампании используют ОДИН И ТОТ ЖЕ base_seed,
%   поэтому реализация i во всех кампаниях получает те же ошибки IMU, ту же
%   выставку и ту же аномалию поля (последняя — из отдельного потока
%   seed_i + GRAV_SEED_OFFSET). Различается только отключённый источник.
%
%   СТАТИСТИКА СЧИТАЕТСЯ ПО ПАРНЫМ РАЗНОСТЯМ:
%       d(i) = E_base(i) - E_without(i)
%   а НЕ по формуле для двух независимых выборок. Разброс выборки при этом
%   сокращается: сравниваются одни и те же реализации, отличающиеся только
%   наличием источника. Это повышает чувствительность в разы при том же N —
%   различия в доли метра становятся разрешимы там, где непарное сравнение
%   потребовало бы сотен реализаций.
%
%   КРИТИЧНО: источник отключается ТОЛЬКО в профиле ИСТИНЫ (prof), а фильтр
%   получает НЕИЗМЕНЁННЫЙ prof_base. Иначе обнуление, например,
%   turnon_bias_sigma обнулило бы и блок P0 для bias — фильтр стал бы
%   абсолютно уверен, что состояния нет, и перестал бы его оценивать.
%   Наблюдался рост КВО на 1.19 м при отключении источника, что абсурдно.
%
%   Вход/выход — см. main_step10_budget.m

    % =====================================================================
    % Шаг Б.1: СПИСОК ИСТОЧНИКОВ
    % =====================================================================
    % Каждый источник задаётся ЯЧЕЙКОЙ путей: один источник может требовать
    % обнуления НЕСКОЛЬКИХ полей.
    %   B2 — уклонение отвеса И аномалия модуля g;
    %   B4 — величина сдвига И флаг включения;
    %   перекос осей — поворот триады (misalign_axes) И неортогональность
    %        (misalign_ortho). Это два разных физических эффекта, но оба
    %        относятся к неидеальности геометрии триады и отключаются вместе;
    %   VRC — достаточно обнулить сам коэффициент, поскольку
    %        b_vrc = VRC * g_rms^2 и при VRC = 0 вклад исчезает независимо
    %        от заданного уровня вибрации.
    %
    % СОСТАВ СВЕРЕН с imu_profile_pulse40_updated:
    % включены все основные моделируемые IMU-источники, B2 и B4.
    % GNSS measurement noise отдельно здесь не изолируется, поскольку
    % sigma_pos/sigma_vel одновременно задают и шум измерения, и R фильтра.
    sources = {
      'turn-on bias акс.',      {'prof.accel.turnon_bias_sigma'}
      'turn-on bias гиро',      {'prof.gyro.turnon_bias_sigma'}
      'in-run bias акс.',       {'prof.accel.inrun_bias_sigma'}
      'in-run bias гиро',       {'prof.gyro.inrun_bias_sigma'}
      'scale factor акс.',      {'prof.accel.scale_factor_sigma'}
      'scale factor гиро',      {'prof.gyro.scale_factor_sigma'}
      'перекос осей акс.',      {'prof.accel.misalign_axes', ...
                                 'prof.accel.misalign_ortho'}
      'перекос осей гиро',      {'prof.gyro.misalign_axes', ...
                                 'prof.gyro.misalign_ortho'}
      'g-sens (ошибка калибр.)',{'prof.gyro.g_sens_cal_sigma'}
      'VRC акс. (вибрация)',    {'prof.accel.vrc'}
      'VRC гиро (вибрация)',    {'prof.gyro.vrc'}
      'шум акс. (VRW)',         {'prof.accel.vrw'}
      'шум гиро (ARW)',         {'prof.gyro.arw'}
      'ошибка выставки',        {'cfg.align_sigma'}
      'B2 аномалия гравитации', {'cfg.defl_vert_sigma', 'cfg.grav_anom_sigma'}
      'B4 сдвиг GNSS',          {'cfg.gnss_time_offset', ...
                                 'cfg.gnss_time_offset_enable'}
    };
    n_src = size(sources,1);

    fprintf('\n=====================================================\n');
    fprintf('  БЮДЖЕТ ОШИБОК: %s\n', prof_base.name);
    fprintf('  %d кампаний x %d реализаций = %d прогонов\n', ...
            n_src+1, N_mc, (n_src+1)*N_mc);
    fprintf('  Сравнение ПАРНОЕ по seed = %d\n', base_seed);
    fprintf('=====================================================\n\n');

    % =====================================================================
    % Шаг Б.2: БАЗОВАЯ КАМПАНИЯ
    % =====================================================================
    fprintf('[база] все источники включены\n');
    mc_base  = run_montecarlo_eskf2(prof_base, prof_base, imu_ideal, truth, c, ...
                                    cfg_base, N_mc, base_seed);
    cep_base = mc_base.cep;
    std_base = mc_base.std_horiz;

    % Базовая выборка горизонтальных ошибок — основа парных разностей
    E_base = mc_base.err_horiz;

    fprintf('  КВО базовый: %.3f м (СКО %.3f м)\n\n', cep_base, std_base);

    % =====================================================================
    % Шаг Б.3: ОТКЛЮЧЕНИЕ ПО ОДНОМУ ИСТОЧНИКУ
    % =====================================================================
    cep_off  = zeros(n_src,1);
    contrib  = zeros(n_src,1);
    signif   = false(n_src,1);
    d_mean   = zeros(n_src,1);   d_std    = zeros(n_src,1);
    d_se     = zeros(n_src,1);   t_stat   = zeros(n_src,1);
    frac_pos = zeros(n_src,1);

    for i = 1:n_src
        name  = sources{i,1};
        paths = sources{i,2};

        % prof — только для ГЕНЕРАЦИИ ошибок; фильтр получает prof_base.
        prof = prof_base;
        cfg  = cfg_base;

        % Обнуляем ВСЕ поля источника. Умножение на ноль сохраняет тип и
        % размерность: для векторов даёт нулевой вектор, для логического
        % флага даёт значение, которое условие трактует как false.
        for ip = 1:numel(paths)
            eval([paths{ip} ' = 0*' paths{ip} ';']);
        end

        fprintf('[%2d/%2d] отключено: %s\n', i, n_src, name);
        mc_i = run_montecarlo_eskf2(prof, prof_base, imu_ideal, truth, c, ...
                                    cfg, N_mc, base_seed);

        cep_off(i) = mc_i.cep;
        contrib(i) = sqrt(max(cep_base^2 - mc_i.cep^2, 0));

        % --- ПАРНАЯ СТАТИСТИКА ---
        d_pair = E_base - mc_i.err_horiz;      % > 0 => источник ухудшал
        okp    = ~isnan(d_pair);
        n_ok   = sum(okp);

        d_mean(i)   = mean(d_pair(okp));
        d_std(i)    = std(d_pair(okp));
        d_se(i)     = d_std(i) / sqrt(max(n_ok,1));
        t_stat(i)   = d_mean(i) / max(d_se(i), eps);
        frac_pos(i) = mean(d_pair(okp) > 0);
        signif(i)   = abs(t_stat(i)) > 2;

        fprintf(['  КВО без него: %.3f м | парная разность %+.4f м ' ...
                 '(СКО ср. %.4f, t = %+.1f) %s\n'], ...
                cep_off(i), d_mean(i), d_se(i), t_stat(i), ...
                ternary_str(signif(i), '[ЗНАЧИМО]', ''));
        fprintf('    ухудшал результат в %.0f%% реализаций\n\n', 100*frac_pos(i));
    end

    % =====================================================================
    % Шаг Б.4: СОРТИРОВКА И УПАКОВКА
    % =====================================================================
    % Сортируем по ПАРНОЙ РАЗНОСТИ: она статистически надёжнее вклада через
    % квадратичное сложение, поскольку не опирается на независимость
    % источников (а они через фильтр коррелированы).
    [~, ord] = sort(d_mean, 'descend');

    budget.names     = sources(ord,1);
    budget.cep_off   = cep_off(ord);
    budget.contrib   = contrib(ord);
    budget.signif    = signif(ord);
    budget.d_mean    = d_mean(ord);
    budget.d_std     = d_std(ord);
    budget.d_se      = d_se(ord);
    budget.t_stat    = t_stat(ord);
    budget.frac_pos  = frac_pos(ord);
    budget.cep_base  = cep_base;
    budget.std_base  = std_base;
    budget.E_base    = E_base;
    budget.mc_base   = mc_base;
    budget.N_mc      = N_mc;
    budget.prof_name = prof_base.name;
    budget.rss       = sqrt(sum(contrib.^2));
    budget.residual  = sqrt(max(cep_base^2 - budget.rss^2, 0));

    % =====================================================================
    % Шаг Б.5: ТАБЛИЦА
    % =====================================================================
    fprintf('\n%s\n', repmat('=',1,94));
    fprintf('  БЮДЖЕТ ОШИБОК: %s (N=%d на кампанию, ПАРНОЕ сравнение)\n', ...
            budget.prof_name, N_mc);
    fprintf('%s\n', repmat('=',1,94));
    fprintf('  %-26s %11s %9s %8s %10s %8s\n', ...
        'Источник', 'dE парн.,м', 'СКО ср.', 't', 'КВО без', 'хуже,%');
    fprintf('  %s\n', repmat('-',1,90));
    for i = 1:n_src
        mark = '  ';
        if budget.signif(i), mark = ' *'; end
        fprintf('  %-26s %11.4f %9.4f %8.1f %10.3f %7.0f%%%s\n', ...
                budget.names{i}, budget.d_mean(i), budget.d_se(i), ...
                budget.t_stat(i), budget.cep_off(i), ...
                100*budget.frac_pos(i), mark);
    end
    fprintf('  %s\n', repmat('-',1,90));
    fprintf('  %-26s %11.3f\n', 'КВО базовый',            cep_base);
    fprintf('%s\n', repmat('=',1,94));
    fprintf('  разность = среднее [E_base(i) - E_without(i)] по парным реализациям\n');
    fprintf('  положительная разность => источник УХУДШАЛ результат\n');
    fprintf('  "хуже,%%" — доля реализаций, где источник ухудшил результат\n');
    fprintf('  * — |t| > 2, разность статистически значима\n');

    n_signif = sum(budget.signif);
    fprintf('\n  Значимых источников: %d из %d\n', n_signif, n_src);
    if n_signif < 3
        % Требуемое N оцениваем по типичному СКО ПАРНОЙ разности: оно много
        % меньше СКО самой выборки, поэтому парный дизайн требует
        % существенно меньше реализаций.
        typ_std = median(d_std(d_std > 0));
        if ~isempty(typ_std) && isfinite(typ_std) && typ_std > 0
            Nreq = ceil((2*typ_std/0.3)^2);
            fprintf('  Для разрешения вкладов уровня 0.3 м нужно N >= %d\n', Nreq);
            fprintf('  (оценка по СКО парной разности %.3f м)\n', typ_std);
        end
    end
end


function s = ternary_str(cond, a, b)
%TERNARY_STR  Выбор строки по условию.
    if cond, s = a; else, s = b; end
end