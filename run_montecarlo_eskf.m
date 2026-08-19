function mc = run_montecarlo_eskf(prof, imu_ideal, truth, c, cfg, N_mc, base_seed)
%RUN_MONTECARLO_ESKF  Монте-Карло полного контура ИНС+GNSS+ESKF.
%
%   Прогоняет N_mc независимых реализаций полного фильтра и собирает
%   статистику ошибки в точке удара. В отличие от run_montecarlo_coast
%   (чистая ИНС без фильтра), здесь работает ВЕСЬ контур, поэтому
%   результат - НАСТОЯЩИЙ КВО с оценкой bias фильтром.
%
%   Что разыгрывается независимо на каждой реализации:
%     - turn-on bias и scale factor IMU (константы на полёт)
%     - in-run bias (процесс Гаусса-Маркова во времени)
%     - шум датчиков (ARW/VRW)
%     - ошибка начальной выставки (крен/тангаж/курс)
%     - шум измерений GNSS
%
%   ВАЖНО о seed: каждая реализация получает СВОЙ seed (base_seed + i),
%   поэтому результаты различаются. Один и тот же base_seed воспроизводит
%   всю кампанию целиком - это нужно для повторяемости результатов.
%
%   Вход:
%     prof       - профиль IMU (например imu_profile_pulse40())
%     imu_ideal  - идеальные показания IMU (из main_step12)
%     truth      - истинная траектория
%     c          - константы
%     cfg        - конфигурация прогона (см. main_step6.m); поле cfg.align_sigma
%                  задаёт СКО ошибки выставки (разыгрывается на каждой реализации)
%     N_mc       - число реализаций
%     base_seed  - базовый seed кампании
%   Выход:
%     mc - структура статистики:
%          .err_horiz  (N_mc x 1) горизонтальная ошибка на ударе [м]
%          .err_3d     (N_mc x 1) полная 3D-ошибка [м]
%          .err_enu    (N_mc x 3) вектор ошибки в ENU [м]
%          .cep        КВО = медиана горизонтальной ошибки [м]
%          .r95        95-й процентиль [м]
%          .ba_resid   (N_mc x 1) остаточная ошибка оценки bias акс. [м/с²]
%          .ba_true    (N_mc x 1) истинный turn-on bias акс. [м/с²]
%          .align_used (N_mc x 3) разыгранная ошибка выставки [рад]
%          .n_fail     число расходившихся прогонов (NaN/Inf)

    % =====================================================================
    % ПРОВЕРКА СОГЛАСОВАННОСТИ КОНФИГУРАЦИИ (частый источник тихих ошибок)
    % =====================================================================
    % Начальная ковариация ориентации P0_att должна быть НЕ МЕНЬШЕ реального
    % СКО ошибки выставки align_sigma. Иначе фильтр "переуверен": он считает
    % ошибку меньше, чем она есть, и недостаточно доверяет измерениям.
    % Требуем запас не менее 1.2x.
    if any(cfg.P0_att(:) < cfg.align_sigma(:))
        warning('run_montecarlo_eskf:P0att', ...
            ['P0_att меньше align_sigma по осям [%s] - фильтр ПЕРЕУВЕРЕН. ' ...
             'Рекомендуется P0_att >= 1.2*align_sigma.'], ...
             num2str(find(cfg.P0_att(:) < cfg.align_sigma(:))'));
    end
    % Проверка достаточности окна GNSS перед коастом: именно там фильтр
    % успевает оценить bias, и от него сильнее всего зависит КВО.
    if size(cfg.outage,1) >= 2
        t_gnss_window = cfg.outage(2,1) - cfg.outage(1,2);
        if t_gnss_window < 10
            warning('run_montecarlo_eskf:shortGNSS', ...
                'Окно GNSS перед коастом всего %.1f с - фильтр может не успеть оценить bias.', ...
                t_gnss_window);
        end
    end

    err_horiz  = zeros(N_mc,1);
    err_3d     = zeros(N_mc,1);
    err_enu    = zeros(N_mc,3);
    ba_resid   = zeros(N_mc,1);
    ba_true    = zeros(N_mc,1);
    align_used = zeros(N_mc,3);
    n_fail     = 0;

    fprintf('Монте-Карло ESKF: %s, N=%d реализаций\n', prof.name, N_mc);
    fprintf('%s\n', repmat('-',1,58));
    t_start = tic;

    for i = 1:N_mc
        seed_i = base_seed + i;

        % --- Розыгрыш ошибок IMU (свои на каждой реализации) ---
        rng_imu  = RandStream('mt19937ar','Seed', seed_i);
        e_imu    = imu_draw_errors(prof, rng_imu);
        imu_meas = add_imu_errors(imu_ideal, prof, e_imu, rng_imu);

        % --- Розыгрыш ошибки начальной выставки ---
        % cfg.align_sigma - СКО [рад] по осям (крен, тангаж, курс)
        cfg_i = cfg;
        cfg_i.align_err = cfg.align_sigma(:) .* randn(rng_imu, 3, 1);
        cfg_i.seed      = seed_i;      % шум GNSS тоже свой

        % --- Прогон полного контура ---
        res = eskf_run(imu_meas, truth, c, prof, cfg_i);

        % --- Сбор результата с проверкой на расходимость ---
        if ~isfinite(res.err_horiz_final) || res.err_horiz_final > 1e5
            n_fail = n_fail + 1;
            err_horiz(i) = NaN;  err_3d(i) = NaN;  err_enu(i,:) = NaN;
        else
            err_horiz(i)  = res.err_horiz_final;
            err_3d(i)     = res.err_3d_final;
            err_enu(i,:)  = res.dr_enu(end,:);
        end
        % Статистика оценки bias: NaN для расходившихся, иначе остаток
        if isnan(err_horiz(i))
            ba_resid(i) = NaN;
        else
            ba_resid(i) = norm(e_imu.accel_turnon_bias - res.nav.ba);
        end
        ba_true(i)      = norm(e_imu.accel_turnon_bias);
        align_used(i,:) = cfg_i.align_err';

        % --- Прогресс каждые 10% ---
        if mod(i, max(1,round(N_mc/10))) == 0 || i == N_mc
            el = toc(t_start);
            fprintf('  %3.0f%%  (%4d/%4d)  прошло %5.1f с, осталось ~%5.1f с\n', ...
                    100*i/N_mc, i, N_mc, el, el*(N_mc-i)/i);
        end
    end

    % =====================================================================
    % Статистика (NaN исключаются)
    % =====================================================================
    ok = ~isnan(err_horiz);
    mc.err_horiz  = err_horiz;
    mc.err_3d     = err_3d;
    mc.err_enu    = err_enu;
    mc.ba_resid   = ba_resid;
    mc.ba_true    = ba_true;
    mc.align_used = align_used;
    mc.n_fail     = n_fail;
    mc.N_mc       = N_mc;
    mc.prof_name  = prof.name;

    mc.cep        = median(err_horiz(ok));      % КВО = медиана
    mc.mean_horiz = mean(err_horiz(ok));
    mc.std_horiz  = std(err_horiz(ok));
    mc.r95        = prctile(err_horiz(ok), 95);
    mc.max_horiz  = max(err_horiz(ok));
    mc.frac_pass  = mean(err_horiz(ok) < cfg.cep_target);

    % Систематическое смещение точки попадания (bias прицеливания)
    mc.bias_enu   = mean(err_enu(ok,:), 1);

    % Доля оценённого bias (только по успешным прогонам)
    ok_b = ~isnan(ba_resid);
    mc.ba_frac_est = 1 - median(ba_resid(ok_b)./ba_true(ok_b));

    fprintf('%s\n', repmat('-',1,58));
    fprintf('Кампания завершена за %.1f с (%.2f с/прогон)\n', toc(t_start), toc(t_start)/N_mc);
end