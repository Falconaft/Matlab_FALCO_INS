function budget = run_error_budget(prof_base, imu_ideal, truth, c, cfg_base, N_mc, base_seed)
%RUN_ERROR_BUDGET  Бюджет ошибок: вклад каждого источника в КВО.
%
%   МЕТОД: базовая кампания (все источники включены), затем по одной кампании
%   на каждый источник с ЭТИМ ИСТОЧНИКОМ ОТКЛЮЧЁННЫМ В ГЕНЕРАЦИИ ОШИБОК.
%
%   КРИТИЧНО (исправлено): источник отключается ТОЛЬКО в профиле ИСТИНЫ.
%   Настройка фильтра (P0, Q) остаётся неизменной. Иначе обнуление, например,
%   turnon_bias_sigma обнуляло бы и блок P0 для bias — фильтр становился бы
%   абсолютно уверен, что bias нет, и переставал его оценивать. Наблюдался
%   рост КВО на 1.19 м при отключении источника, что физически невозможно.
%
%   Вклад источника:
%       вклад_i = sqrt( max(КВО_база² - КВО_без_i², 0) )
%   Корректно при независимых источниках, складывающихся в квадратуре.
%
%   СТАТИСТИЧЕСКАЯ ЗНАЧИМОСТЬ: разрешаемый вклад ограничен размером выборки.
%   Порог различимости двух кампаний ~2*sqrt(2)*1.253*sigma/sqrt(N).
%   При sigma~1.1 м: N=10 разрешает только вклады >1.2 м, N=60 -> >0.5 м,
%   N=170 -> >0.3 м. Скрипт печатает порог и помечает незначимые вклады.
%
%   Вход/выход — см. main_step10_budget.m

    % =====================================================================
    % Шаг Б.1: СПИСОК ИСТОЧНИКОВ
    % =====================================================================
    sources = {
      'turn-on bias акс.',      'prof.accel.turnon_bias_sigma'
      'turn-on bias гиро',      'prof.gyro.turnon_bias_sigma'
      'in-run bias акс.',       'prof.accel.inrun_bias_sigma'
      'in-run bias гиро',       'prof.gyro.inrun_bias_sigma'
      'scale factor акс.',      'prof.accel.scale_factor_sigma'
      'scale factor гиро',      'prof.gyro.scale_factor_sigma'
      'g-sens (ошибка калибр.)','prof.gyro.g_sens_cal_sigma'
      'шум акс. (VRW)',         'prof.accel.vrw'
      'шум гиро (ARW)',         'prof.gyro.arw'
      'ошибка выставки',        'cfg.align_sigma'
      'аномалия гравитации',    'cfg.defl_vert_sigma'
    };
    n_src = size(sources,1);

    fprintf('\n=====================================================\n');
    fprintf('  БЮДЖЕТ ОШИБОК: %s\n', prof_base.name);
    fprintf('  %d кампаний x %d реализаций = %d прогонов\n', ...
            n_src+1, N_mc, (n_src+1)*N_mc);
    fprintf('=====================================================\n\n');

    % =====================================================================
    % Шаг Б.2: БАЗОВАЯ КАМПАНИЯ
    % =====================================================================
    fprintf('[база] все источники включены\n');
    mc_base  = run_montecarlo_eskf2(prof_base, prof_base, imu_ideal, truth, c, ...
                                    cfg_base, N_mc, base_seed);
    cep_base = mc_base.cep;
    std_base = mc_base.std_horiz;

    % Порог статистической различимости (2 sigma на разность двух медиан)
    se_med  = 1.253*std_base/sqrt(N_mc);
    thr_sig = 2*sqrt(2)*se_med;

    fprintf('  КВО базовый: %.3f м (СКО %.3f м)\n', cep_base, std_base);
    fprintf('  Порог значимости разности: %.3f м\n\n', thr_sig);

    % =====================================================================
    % Шаг Б.3: ОТКЛЮЧЕНИЕ ПО ОДНОМУ ИСТОЧНИКУ
    % =====================================================================
    cep_off = zeros(n_src,1);
    contrib = zeros(n_src,1);
    drop    = zeros(n_src,1);
    signif  = false(n_src,1);

    for i = 1:n_src
        name = sources{i,1};
        path = sources{i,2};

        % ВАЖНО: prof — только для ГЕНЕРАЦИИ ошибок (истина).
        % prof_base передаётся фильтру НЕИЗМЕННЫМ.
        prof = prof_base;
        cfg  = cfg_base;
        eval([path ' = 0*' path ';']);

        fprintf('[%2d/%2d] отключено: %s\n', i, n_src, name);
        mc_i = run_montecarlo_eskf2(prof, prof_base, imu_ideal, truth, c, ...
                                    cfg, N_mc, base_seed);

        cep_off(i) = mc_i.cep;
        drop(i)    = cep_base - mc_i.cep;
        contrib(i) = sqrt(max(cep_base^2 - mc_i.cep^2, 0));
        signif(i)  = abs(drop(i)) > thr_sig;

        fprintf('  КВО без него: %.3f м  (снижение %+.3f м, вклад %.3f м) %s\n\n', ...
                cep_off(i), drop(i), contrib(i), ...
                char("" + string(repmat('*',1,double(signif(i))))));
    end

    % =====================================================================
    % Шаг Б.4: СОРТИРОВКА И УПАКОВКА
    % =====================================================================
    [~, ord] = sort(contrib, 'descend');
    budget.names     = sources(ord,1);
    budget.cep_off   = cep_off(ord);
    budget.contrib   = contrib(ord);
    budget.drop      = drop(ord);
    budget.signif    = signif(ord);
    budget.cep_base  = cep_base;
    budget.std_base  = std_base;
    budget.thr_sig   = thr_sig;
    budget.mc_base   = mc_base;
    budget.N_mc      = N_mc;
    budget.prof_name = prof_base.name;
    budget.rss       = sqrt(sum(contrib.^2));
    budget.residual  = sqrt(max(cep_base^2 - budget.rss^2, 0));

    % =====================================================================
    % Шаг Б.5: ТАБЛИЦА
    % =====================================================================
    fprintf('\n==================================================================\n');
    fprintf('  БЮДЖЕТ ОШИБОК: %s (N=%d на кампанию)\n', budget.prof_name, N_mc);
    fprintf('==================================================================\n');
    fprintf('  %-26s %9s %9s %8s %6s\n', 'Источник','вклад,м','КВО без','доля','знач.');
    fprintf('  %s\n', repmat('-',1,62));
    for i = 1:n_src
        frac = 100*budget.contrib(i)^2 / max(cep_base^2, eps);
        mark = '   ';
        if budget.signif(i), mark = ' * '; end
        fprintf('  %-26s %9.3f %9.3f %7.1f%% %5s\n', ...
                budget.names{i}, budget.contrib(i), budget.cep_off(i), frac, mark);
    end
    fprintf('  %s\n', repmat('-',1,62));
    fprintf('  %-26s %9.3f\n', 'КВО базовый', cep_base);
    fprintf('  %-26s %9.3f\n', 'корень суммы квадратов', budget.rss);
    fprintf('  %-26s %9.3f\n', 'необъяснённый остаток', budget.residual);
    fprintf('  %-26s %9.3f\n', 'порог значимости', thr_sig);
    fprintf('==================================================================\n');
    fprintf('  * — вклад статистически значим при данном N\n');

    n_signif = sum(budget.signif);
    if n_signif < 3
        Nreq = ceil((2*sqrt(2)*1.253*std_base/0.5)^2);
        fprintf('\n  ВНИМАНИЕ: значимы только %d источник(ов) из %d.\n', n_signif, n_src);
        fprintf('  Для разрешения вкладов уровня 0.5 м нужно N >= %d.\n', Nreq);
    end
end