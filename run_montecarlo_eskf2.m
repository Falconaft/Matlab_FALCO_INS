function mc = run_montecarlo_eskf2(prof_truth, prof_filter, imu_ideal, truth, c, cfg, N_mc, base_seed)
%RUN_MONTECARLO_ESKF2  Монте-Карло с РАЗДЕЛЬНЫМИ профилями истины и фильтра.
%
%   ПАРАЛЛЕЛЬНАЯ ВЕРСИЯ ДЛЯ PARFOR.
%
%   Основа: актуальный run_montecarlo_eskf2.m из ветки
%   feature/consistency-diagnostics.
%
%   Главное отличие от последовательной версии:
%     - независимые реализации Monte-Carlo выполняются через parfor;
%     - диагностические временные ряды NEES/NIS/rho внутри parfor
%       собираются через cell-массивы;
%     - после parfor они упаковываются в обычные матрицы;
%     - seed каждой реализации НЕ изменён:
%           seed_i = base_seed + i
%       поэтому парность между кампаниями бюджета ошибок сохраняется.
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
%     mc - статистика (та же структура, что у последовательной версии)

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

    % КАЧЕСТВО ОЦЕНКИ BIAS — строго в момент ПЕРЕД терминальным коастом.
    ba_true_pre  = nan(N_mc,3);   ba_est_pre  = nan(N_mc,3);
    bg_true_pre  = nan(N_mc,3);   bg_est_pre  = nan(N_mc,3);
    ba_resid_pre = nan(N_mc,1);   bg_resid_pre = nan(N_mc,1);

    % =====================================================================
    % МЕТРИКИ СОСТОЯТЕЛЬНОСТИ ФИЛЬТРА (NEES / NIS)
    % =====================================================================
    % В последовательной версии матрицы NEES_mat/NIS_mat создавались
    % динамически внутри первой итерации:
    %
    %   if isempty(NEES_mat)
    %       NEES_mat = nan(...);
    %   end
    %
    % Такая схема несовместима с parfor: worker не может безопасно менять
    % размер общей переменной в зависимости от того, какая итерация
    % завершилась первой.
    %
    % Поэтому каждый worker сохраняет свой временной ряд в ОТДЕЛЬНУЮ
    % ячейку. После завершения parfor ряды собираются в матрицы.
    NEES_cell  = cell(N_mc,1);
    NEESR_cell = cell(N_mc,1);
    NEESV_cell = cell(N_mc,1);
    NIS_cell   = cell(N_mc,1);
    RHOP_cell  = cell(N_mc,1);
    RHOV_cell  = cell(N_mc,1);
    TDIAG_cell = cell(N_mc,1);
    TGNSS_cell = cell(N_mc,1);

    % Итоговые матрицы будут сформированы ПОСЛЕ parfor.
    NEES_mat  = [];
    NEESR_mat = [];
    NEESV_mat = [];
    NIS_mat   = [];
    RHOP_mat  = [];
    RHOV_mat  = [];
    t_diag    = [];
    t_gnss    = [];

    % Усреднённые по времени величины — ТОЛЬКО как компактная сводка.
    % Строгим 95% тестом они НЕ являются.
    nees_gnss  = nan(N_mc,1);
    nees_coast = nan(N_mc,1);
    nis_tavg   = nan(N_mc,1);

    % МАСШТАБНЫЕ КОЭФФИЦИЕНТЫ: истина и оценка в ДВУХ точках.
    bg_true_pre_ax = nan(N_mc,3);   bg_est_pre_ax = nan(N_mc,3);
    bg_true_end_ax = nan(N_mc,3);   bg_est_end_ax = nan(N_mc,3);
    rho_pre_ax     = nan(N_mc,3);

    sg_true_all = nan(N_mc,3);   sa_true_all = nan(N_mc,3);
    sg_pre      = nan(N_mc,3);   sa_pre      = nan(N_mc,3);
    sg_end      = nan(N_mc,3);   sa_end      = nan(N_mc,3);

    % Состояние непосредственно перед терминальным коастом
    pre_horiz  = nan(N_mc,1);
    pre_dv     = nan(N_mc,1);
    pre_dpsi   = nan(N_mc,1);
    align_used = zeros(N_mc,3);

    % Применённые начальные ошибки (диагностика)
    init_pos_used = zeros(N_mc,3);
    init_vel_used = zeros(N_mc,3);

    % Для parfor вместо reduction-переменной n_fail используем sliced-флаг.
    fail_flag = false(N_mc,1);

    % Диагностика фактически разыгранных аномалий поля
    dg_info = struct('defl_E_arcsec',cell(N_mc,1), 'defl_N_arcsec',[], ...
                     'defl_total_arcsec',[], 'anom_vert_mgal',[], ...
                     'dg_horiz_mps2',[], 'dg_norm_mps2',[]);
    dg_vec  = zeros(N_mc,3);

    % Начало терминального коаста
    if size(cfg.outage,1) >= 2
        t_coast_start = cfg.outage(2,1);
    else
        t_coast_start = truth.t(end);
    end

    % =====================================================================
    % ПАРАЛЛЕЛЬНЫЙ MONTE-CARLO
    % =====================================================================
    % MATLAB автоматически поднимет Parallel Pool при первом parfor, если
    % пул ещё не запущен. Порядок выполнения worker-ами не влияет на
    % воспроизводимость, потому что seed задаётся явно через индекс i.
    parfor i = 1:N_mc
        seed_i = base_seed + i;

        % Локальные диагностические ряды ЭТОЙ реализации.
        % Они присваиваются cell-массивам только один раз в конце итерации.
        nees_i   = [];
        neesr_i  = [];
        neesv_i  = [];
        nis_i    = [];
        rhop_i   = [];
        rhov_i   = [];
        tdiag_i  = [];
        tgnss_i  = [];

        % --- Розыгрыш ошибок по профилю ИСТИНЫ ---
        rng_imu  = RandStream('mt19937ar','Seed', seed_i);
        e_imu    = imu_draw_errors(prof_truth, rng_imu);
        imu_meas = add_imu_errors(imu_ideal, prof_truth, e_imu, rng_imu);

        % --- Ошибка выставки (своя на каждой реализации) ---
        cfg_i = cfg;
        cfg_i.align_err = cfg.align_sigma(:) .* randn(rng_imu, 3, 1);

        % НАЧАЛЬНАЯ ОШИБКА ПОЗИЦИИ И СКОРОСТИ
        if isfield(cfg,'init_pos_sigma')
            cfg_i.init_pos_err = cfg.init_pos_sigma(:) .* randn(rng_imu, 3, 1);
        else
            cfg_i.init_pos_err = zeros(3,1);
        end

        if isfield(cfg,'init_vel_sigma')
            cfg_i.init_vel_err = cfg.init_vel_sigma(:) .* randn(rng_imu, 3, 1);
        else
            cfg_i.init_vel_err = zeros(3,1);
        end

        cfg_i.seed = seed_i;

        % Аномалия гравитации: отдельный поток.
        rng_grav = RandStream('mt19937ar','Seed', seed_i + GRAV_SEED_OFFSET);
        [cfg_i.dg_model, info_i] = draw_gravity_anomaly(cfg, truth.Cen, rng_grav);

        dg_info(i)  = info_i;
        dg_vec(i,:) = cfg_i.dg_model';

        % --- Прогон фильтра, настроенного по профилю ФИЛЬТРА ---
        res = eskf_run(imu_meas, truth, c, prof_filter, cfg_i);

        % --- Сбор результата ---
        if ~isfinite(res.err_horiz_final) || res.err_horiz_final > 1e5
            fail_flag(i) = true;

            err_horiz(i) = NaN;
            err_3d(i)    = NaN;
            err_enu(i,:) = NaN;

            ba_resid_pre(i) = NaN;
            bg_resid_pre(i) = NaN;
        else
            err_horiz(i) = res.err_horiz_final;
            err_3d(i)    = res.err_3d_final;
            err_enu(i,:) = res.dr_enu(end,:);

            % Последняя диагностическая точка ДО начала аутэйджа
            k_pre = find(res.t < t_coast_start, 1, 'last');
            if isempty(k_pre)
                k_pre = 1;
            end

            % =============================================================
            % Метрики состоятельности
            % =============================================================
            if isfield(res,'nees6')
                tdiag_i = res.t(:)';
                nees_i  = res.nees6(:)';

                if isfield(res,'nees_r')
                    neesr_i = res.nees_r(:)';
                    neesv_i = res.nees_v(:)';
                end

                % Компактная сводка по участкам (НЕ строгий тест)
                i_coast = res.t >= t_coast_start;
                i_gnss  = ~i_coast & (res.t > cfg.outage(1,2));

                if any(i_gnss)
                    nees_gnss(i) = mean(res.nees6(i_gnss), 'omitnan');
                end
                if any(i_coast)
                    nees_coast(i) = mean(res.nees6(i_coast), 'omitnan');
                end
            end

            if isfield(res,'gnss_nis') && ~isempty(res.gnss_nis)
                tgnss_i = res.gnss_t(:)';
                nis_i   = res.gnss_nis(:)';

                if isfield(res,'gnss_rho_pos')
                    rhop_i = mean(res.gnss_rho_pos, 2)';
                    rhov_i = mean(res.gnss_rho_vel, 2)';
                end

                nis_tavg(i) = mean(res.gnss_nis, 'omitnan');
            end

            % --- Смещение гироскопа по осям в двух точках ---
            if isfield(res,'bg_true')
                bg_true_pre_ax(i,:) = res.bg_true(k_pre,:);
                bg_est_pre_ax(i,:)  = res.bg_est(k_pre,:);
                bg_true_end_ax(i,:) = res.bg_true(end,:);
                bg_est_end_ax(i,:)  = res.bg_est(end,:);
            end

            if isfield(res,'rho_bg_sg')
                rho_pre_ax(i,:) = res.rho_bg_sg(k_pre,:);
            end

            % --- Масштабные коэффициенты в двух точках ---
            if isfield(res,'sg_est')
                sg_true_all(i,:) = e_imu.gyro_scale(:)';
                sa_true_all(i,:) = e_imu.accel_scale(:)';
                sg_pre(i,:)      = res.sg_est(k_pre,:);
                sa_pre(i,:)      = res.sa_est(k_pre,:);
                sg_end(i,:)      = res.sg_est(end,:);
                sa_end(i,:)      = res.sa_est(end,:);
            end

            pre_horiz(i) = hypot(res.dr_ned(k_pre,1), res.dr_ned(k_pre,2));
            pre_dv(i)    = norm(res.dv(k_pre,:));
            pre_dpsi(i)  = norm(res.dpsi(k_pre,:)) * 1e3;   % мрад

            % Истинное ПОЛНОЕ смещение в тот же момент
            [~, k_truth] = min(abs(imu_meas.t - res.t(k_pre)));

            bg_true_pre(i,:) = imu_meas.bg_total(k_truth,:);
            ba_true_pre(i,:) = imu_meas.ba_total(k_truth,:);
            bg_est_pre(i,:)  = res.bg_est(k_pre,:);
            ba_est_pre(i,:)  = res.ba_est(k_pre,:);

            bg_resid_pre(i) = norm(bg_true_pre(i,:) - bg_est_pre(i,:));
            ba_resid_pre(i) = norm(ba_true_pre(i,:) - ba_est_pre(i,:));
        end

        align_used(i,:)    = cfg_i.align_err';
        init_pos_used(i,:) = cfg_i.init_pos_err';
        init_vel_used(i,:) = cfg_i.init_vel_err';

        % Единственное sliced-присваивание диагностических cell-массивов.
        NEES_cell{i}  = nees_i;
        NEESR_cell{i} = neesr_i;
        NEESV_cell{i} = neesv_i;
        NIS_cell{i}   = nis_i;
        RHOP_cell{i}  = rhop_i;
        RHOV_cell{i}  = rhov_i;
        TDIAG_cell{i} = tdiag_i;
        TGNSS_cell{i} = tgnss_i;
    end

    n_fail = sum(fail_flag);

    % =====================================================================
    % СБОР ДИАГНОСТИЧЕСКИХ РЯДОВ ПОСЛЕ PARFOR
    % =====================================================================
    % Здесь мы уже снова в клиентском MATLAB, поэтому можно безопасно
    % определить размер временной сетки по первой успешной реализации.
    [NEES_mat, t_diag] = pack_series_with_time(NEES_cell, TDIAG_cell, N_mc);
    NEESR_mat = pack_series(NEESR_cell, N_mc, size(NEES_mat,2));
    NEESV_mat = pack_series(NEESV_cell, N_mc, size(NEES_mat,2));

    [NIS_mat, t_gnss] = pack_series_with_time(NIS_cell, TGNSS_cell, N_mc);
    RHOP_mat = pack_series(RHOP_cell, N_mc, size(NIS_mat,2));
    RHOV_mat = pack_series(RHOV_cell, N_mc, size(NIS_mat,2));

    % =====================================================================
    % СТАТИСТИКА
    % =====================================================================
    ok = ~isnan(err_horiz);

    mc.err_horiz = err_horiz;
    mc.err_3d    = err_3d;
    mc.err_enu   = err_enu;

    mc.align_used    = align_used;
    mc.init_pos_used = init_pos_used;
    mc.init_vel_used = init_vel_used;

    % Качество оценки bias перед коастом
    mc.ba_true_pre  = ba_true_pre;
    mc.ba_est_pre   = ba_est_pre;
    mc.bg_true_pre  = bg_true_pre;
    mc.bg_est_pre   = bg_est_pre;
    mc.ba_resid_pre = ba_resid_pre;
    mc.bg_resid_pre = bg_resid_pre;

    mc.pre_horiz = pre_horiz;
    mc.pre_dv    = pre_dv;
    mc.pre_dpsi  = pre_dpsi;

    mc.t_coast_start = t_coast_start;
    mc.n_fail        = n_fail;
    mc.N_mc          = N_mc;
    mc.prof_name     = prof_truth.name;

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

    mc.bias_enu = mean(err_enu(ok,:), 1);

    % --- Диагностика аномалии гравитации ---
    mc.dg_info = dg_info;
    mc.dg_vec  = dg_vec;
    mc.grav_seed_offset = GRAV_SEED_OFFSET;

    % =====================================================================
    % СТАТИСТИКА СОСТОЯТЕЛЬНОСТИ: ANEES / ANIS
    % =====================================================================
    dim = 6;
    mc.consist_dim = dim;

    if ~isempty(NEES_mat)
        n_ok_k = sum(~isnan(NEES_mat), 1);

        mc.anees   = mean(NEES_mat, 1, 'omitnan');
        mc.anees_t = t_diag;
        mc.anees_n = n_ok_k;

        mc.anees_lo = nan(size(mc.anees));
        mc.anees_hi = nan(size(mc.anees));

        for k = 1:numel(mc.anees)
            nk = n_ok_k(k);
            if nk > 0
                [lo, hi] = chi2_bounds_wh(nk*dim, 0.95);
                mc.anees_lo(k) = lo/nk;
                mc.anees_hi(k) = hi/nk;
            end
        end

        mc.NEES_mat = NEES_mat;

        if any(n_ok_k > 0)
            mc.consist_n_typ = median(n_ok_k(n_ok_k > 0));
        else
            mc.consist_n_typ = NaN;
        end
    end

    % --- ANEES по каналам, dim = 3 ---
    dim3 = 3;

    if ~isempty(NEESR_mat)
        for tag = {'r','v'}
            t_ = tag{1};

            if strcmp(t_,'r')
                M = NEESR_mat;
            else
                M = NEESV_mat;
            end

            n_ok_c = sum(~isnan(M), 1);
            a      = mean(M, 1, 'omitnan');

            lo_ = nan(size(a));
            hi_ = nan(size(a));

            for k = 1:numel(a)
                nk = n_ok_c(k);
                if nk > 0
                    [l, h] = chi2_bounds_wh(nk*dim3, 0.95);
                    lo_(k) = l/nk;
                    hi_(k) = h/nk;
                end
            end

            mc.(['anees_' t_])       = a;
            mc.(['anees_' t_ '_lo']) = lo_;
            mc.(['anees_' t_ '_hi']) = hi_;
        end

        mc.anees_rv_dim = dim3;
        mc.NEESR_mat    = NEESR_mat;
        mc.NEESV_mat    = NEESV_mat;
    end

    if ~isempty(NIS_mat)
        n_ok_g = sum(~isnan(NIS_mat), 1);

        mc.anis   = mean(NIS_mat, 1, 'omitnan');
        mc.anis_t = t_gnss;
        mc.anis_n = n_ok_g;

        mc.anis_lo = nan(size(mc.anis));
        mc.anis_hi = nan(size(mc.anis));

        for k = 1:numel(mc.anis)
            nk = n_ok_g(k);
            if nk > 0
                [lo, hi] = chi2_bounds_wh(nk*dim, 0.95);
                mc.anis_lo(k) = lo/nk;
                mc.anis_hi(k) = hi/nk;
            end
        end

        mc.NIS_mat = NIS_mat;
    end

    % --- Доминирование R в инновационной ковариации, ПО КАНАЛАМ ---
    if ~isempty(RHOP_mat)
        mc.rho_pos = mean(RHOP_mat, 1, 'omitnan');
        mc.rho_vel = mean(RHOV_mat, 1, 'omitnan');
        mc.rho_t   = t_gnss;

        finite_pos = mc.rho_pos(isfinite(mc.rho_pos));
        finite_vel = mc.rho_vel(isfinite(mc.rho_vel));

        if isempty(finite_pos)
            mc.rho_pos_med = NaN;
        else
            mc.rho_pos_med = median(finite_pos);
        end

        if isempty(finite_vel)
            mc.rho_vel_med = NaN;
        else
            mc.rho_vel_med = median(finite_vel);
        end

        mc.RHOP_mat = RHOP_mat;
        mc.RHOV_mat = RHOV_mat;
    end

    % Компактные средние по времени
    mc.nees_gnss  = nees_gnss;
    mc.nees_coast = nees_coast;
    mc.nis_tavg   = nis_tavg;

    mc.nees_gnss_mean  = mean(nees_gnss(~isnan(nees_gnss)));
    mc.nees_coast_mean = mean(nees_coast(~isnan(nees_coast)));
    mc.nis_tavg_mean   = mean(nis_tavg(~isnan(nis_tavg)));

    % --- Метрики оценки смещения гироскопа (в °/ч) ---
    DPH = 180/pi*3600;

    if any(isfinite(bg_true_pre_ax(:)))
        mc.gb.pre = estimation_quality_metrics( ...
            bg_true_pre_ax, bg_est_pre_ax, DPH);

        mc.gb.end = estimation_quality_metrics( ...
            bg_true_end_ax, bg_est_end_ax, DPH);

        mc.gb.rho_pre_med = median(rho_pre_ax, 1, 'omitnan');

        mc.gb.bg_true_pre = bg_true_pre_ax;
        mc.gb.bg_est_pre  = bg_est_pre_ax;
        mc.gb.bg_true_end = bg_true_end_ax;
        mc.gb.bg_est_end  = bg_est_end_ax;
        mc.gb.rho_pre     = rho_pre_ax;
    end

    % --- Метрики оценки масштабных коэффициентов ---
    if any(isfinite(sg_true_all(:)))
        mc.sf.gyro_pre  = scale_factor_metrics(sg_true_all, sg_pre);
        mc.sf.gyro_end  = scale_factor_metrics(sg_true_all, sg_end);
        mc.sf.accel_pre = scale_factor_metrics(sa_true_all, sa_pre);
        mc.sf.accel_end = scale_factor_metrics(sa_true_all, sa_end);

        mc.sf.sg_true = sg_true_all;
        mc.sf.sa_true = sa_true_all;
        mc.sf.sg_pre  = sg_pre;
        mc.sf.sa_pre  = sa_pre;
        mc.sf.sg_end  = sg_end;
        mc.sf.sa_end  = sa_end;
    end

    mc.pre_horiz_med = median(pre_horiz(~isnan(pre_horiz)));
    mc.pre_dv_med    = median(pre_dv(~isnan(pre_dv)));
    mc.pre_dpsi_med  = median(pre_dpsi(~isnan(pre_dpsi)));

    % Медианы истинного смещения и остатка оценки перед коастом
    mc.bg_true_pre_med  = median(vecnorm( ...
        bg_true_pre(~isnan(bg_resid_pre),:), 2, 2));

    mc.ba_true_pre_med  = median(vecnorm( ...
        ba_true_pre(~isnan(ba_resid_pre),:), 2, 2));

    mc.bg_resid_pre_med = median(bg_resid_pre(~isnan(bg_resid_pre)));
    mc.ba_resid_pre_med = median(ba_resid_pre(~isnan(ba_resid_pre)));
end


function [M, t] = pack_series_with_time(series_cell, time_cell, N_rows)
%PACK_SERIES_WITH_TIME
%   Собирает временные ряды отдельных реализаций после parfor.
%   Первая непустая реализация определяет ожидаемую длину ряда.
%
%   Если у отдельной реализации длина отличается, её строка остаётся NaN.
%   Это повторяет смысл проверки размеров из исходной последовательной
%   версии, но без динамического изменения общей матрицы внутри parfor.

    M = [];
    t = [];

    first = find(~cellfun(@isempty, series_cell), 1, 'first');
    if isempty(first)
        return;
    end

    n = numel(series_cell{first});
    M = nan(N_rows, n);

    if ~isempty(time_cell{first})
        t = time_cell{first}(:)';
    end

    for i = 1:N_rows
        row = series_cell{i};

        if ~isempty(row) && numel(row) == n
            M(i,:) = row(:)';
        end
    end
end


function M = pack_series(series_cell, N_rows, expected_n)
%PACK_SERIES
%   Собирает ряды без отдельной временной сетки.
%
%   expected_n обычно берётся из уже собранной NEES_mat или NIS_mat.

    M = [];

    if expected_n <= 0
        return;
    end

    if ~any(~cellfun(@isempty, series_cell))
        return;
    end

    M = nan(N_rows, expected_n);

    for i = 1:N_rows
        row = series_cell{i};

        if ~isempty(row) && numel(row) == expected_n
            M(i,:) = row(:)';
        end
    end
end


function [lo, hi] = chi2_bounds_wh(k, conf)
%CHI2_BOUNDS_WH  Квантили распределения хи-квадрат без Statistics Toolbox.
%
%   Аппроксимация Уилсона-Хилферти:
%       chi2_p(k) ~ k * (1 - 2/(9k) + z_p*sqrt(2/(9k)))^3
%   где z_p — квантиль стандартного нормального распределения.
%
%   Точность при k >= 30 лучше 0.1%, что для наших k = N*6 >= 60 избыточно.
%   Использована вместо chi2inv, поскольку Statistics Toolbox может
%   отсутствовать.
%
%   Вход:  k    - число степеней свободы
%          conf - доверительная вероятность (например 0.95)
%   Выход: lo, hi - нижняя и верхняя границы

    alpha = 1 - conf;

    % Квантили нормального распределения для двустороннего интервала.
    z = sqrt(2) * erfinv(1 - alpha);

    lo = k * (1 - 2/(9*k) - z*sqrt(2/(9*k)))^3;
    hi = k * (1 - 2/(9*k) + z*sqrt(2/(9*k)))^3;
end
