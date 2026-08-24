function mc = run_montecarlo_eskf2(prof_truth, prof_filter, imu_ideal, truth, c, cfg, N_mc, base_seed)
%RUN_MONTECARLO_ESKF2  Монте-Карло с РАЗДЕЛЬНЫМИ профилями истины и фильтра.
%
%   Отличие от run_montecarlo_eskf: профиль, по которому ГЕНЕРИРУЮТСЯ ошибки
%   датчика (prof_truth), отделён от профиля, по которому НАСТРАИВАЕТСЯ фильтр
%   (prof_filter: начальная ковариация P0 и матрица шума Q).
%
%   Зачем это нужно для бюджета ошибок:
%     Обнуление параметра в едином профиле отключало источник ошибки И
%     одновременно обнуляло соответствующий блок P0 — фильтр становился
%     абсолютно уверен, что этого состояния нет, и переставал его оценивать.
%     Из-за этого отключение turn-on bias гироскопа УХУДШАЛО КВО на 1.19 м,
%     что физически невозможно.
%
%     Правильный эксперимент: убрать источник из РЕАЛЬНОСТИ, оставив фильтр
%     настроенным как обычно (он не знает, что источника нет).
%
%   Также позволяет исследовать РАССОГЛАСОВАНИЕ настройки: например, задать
%   фильтру завышенную/заниженную оценку шумов и посмотреть на робастность.
%
%   Вход:
%     prof_truth  - профиль для генерации ошибок датчика (что происходит в реальности)
%     prof_filter - профиль для настройки фильтра (P0, Q) — что фильтр предполагает
%     imu_ideal   - идеальные показания IMU
%     truth       - траектория
%     c           - константы
%     cfg         - конфигурация прогона (cfg.align_sigma — СКО ошибки выставки)
%     N_mc        - число реализаций
%     base_seed   - базовый seed кампании
%   Выход:
%     mc - статистика (та же структура, что у run_montecarlo_eskf)

    % =====================================================================
    % ОТДЕЛЬНЫЙ ПОТОК ДЛЯ АНОМАЛИИ ГРАВИТАЦИИ
    % =====================================================================
    % Аномалия разыгрывается из СВОЕГО потока, а не из того, где берутся
    % ошибки IMU и выставка. Зачем:
    %   - реализация поля воспроизводима независимо от модели датчика;
    %   - изменение числа обращений к randn в imu_draw_errors/add_imu_errors
    %     НЕ сдвигает разыгранную аномалию;
    %   - B2 можно независимо включать и выключать в бюджете ошибок,
    %     не меняя прочие реализации.
    % Смещение выбрано заведомо большим числа реализаций любой кампании,
    % чтобы потоки seed_i и seed_i+OFFSET не пересеклись.
    GRAV_SEED_OFFSET = 100000;

    err_horiz  = zeros(N_mc,1);
    err_3d     = zeros(N_mc,1);
    err_enu    = zeros(N_mc,3);
    ba_resid   = zeros(N_mc,1);
    ba_true    = zeros(N_mc,1);
    align_used = zeros(N_mc,3);
    n_fail     = 0;

    % Диагностика фактически разыгранных аномалий поля
    dg_info = struct('defl_E_arcsec',cell(N_mc,1), 'defl_N_arcsec',[], ...
                     'defl_total_arcsec',[], 'anom_vert_mgal',[], ...
                     'dg_horiz_mps2',[], 'dg_norm_mps2',[]);
    dg_vec  = zeros(N_mc,3);

    for i = 1:N_mc
        seed_i = base_seed + i;

        % --- Розыгрыш ошибок по профилю ИСТИНЫ ---
        rng_imu  = RandStream('mt19937ar','Seed', seed_i);
        e_imu    = imu_draw_errors(prof_truth, rng_imu);
        imu_meas = add_imu_errors(imu_ideal, prof_truth, e_imu, rng_imu);

        % --- Ошибка выставки (своя на каждой реализации) ---
        cfg_i = cfg;
        cfg_i.align_err = cfg.align_sigma(:) .* randn(rng_imu, 3, 1);
        cfg_i.seed      = seed_i;

        % Аномалия гравитации: своя на каждой реализации (иначе она стала бы
        % детерминированной и дала бы ложную систематику по ансамблю) и из
        % ОТДЕЛЬНОГО потока (см. пояснение к GRAV_SEED_OFFSET выше).
        rng_grav = RandStream('mt19937ar','Seed', seed_i + GRAV_SEED_OFFSET);
        [cfg_i.dg_model, info_i] = draw_gravity_anomaly(cfg, truth.Cen, rng_grav);
        dg_info(i) = info_i;
        dg_vec(i,:) = cfg_i.dg_model';

        % --- Прогон фильтра, настроенного по профилю ФИЛЬТРА ---
        res = eskf_run(imu_meas, truth, c, prof_filter, cfg_i);

        % --- Сбор результата ---
        if ~isfinite(res.err_horiz_final) || res.err_horiz_final > 1e5
            n_fail = n_fail + 1;
            err_horiz(i) = NaN;  err_3d(i) = NaN;  err_enu(i,:) = NaN;
            ba_resid(i)  = NaN;
        else
            err_horiz(i) = res.err_horiz_final;
            err_3d(i)    = res.err_3d_final;
            err_enu(i,:) = res.dr_enu(end,:);
            ba_resid(i)  = norm(e_imu.accel_turnon_bias - res.nav.ba);
        end
        ba_true(i)      = norm(e_imu.accel_turnon_bias);
        align_used(i,:) = cfg_i.align_err';
    end

    % --- Статистика ---
    ok = ~isnan(err_horiz);
    mc.err_horiz  = err_horiz;    mc.err_3d  = err_3d;    mc.err_enu = err_enu;
    mc.ba_resid   = ba_resid;     mc.ba_true = ba_true;   mc.align_used = align_used;
    mc.n_fail     = n_fail;       mc.N_mc    = N_mc;
    mc.prof_name  = prof_truth.name;

    mc.cep        = median(err_horiz(ok));
    mc.mean_horiz = mean(err_horiz(ok));
    mc.std_horiz  = std(err_horiz(ok));
    mc.r95        = prctile(err_horiz(ok), 95);
    mc.max_horiz  = max(err_horiz(ok));
    if isfield(cfg,'cep_target')
        mc.frac_pass = mean(err_horiz(ok) < cfg.cep_target);
    else
        mc.frac_pass = NaN;
    end
    mc.bias_enu   = mean(err_enu(ok,:), 1);

    % --- Диагностика аномалии гравитации ---
    mc.dg_info = dg_info;
    mc.dg_vec  = dg_vec;
    mc.grav_seed_offset = GRAV_SEED_OFFSET;

    ok_b = ~isnan(ba_resid) & (ba_true > 0);
    if any(ok_b)
        mc.ba_frac_est = 1 - median(ba_resid(ok_b)./ba_true(ok_b));
    else
        mc.ba_frac_est = NaN;
    end
end
