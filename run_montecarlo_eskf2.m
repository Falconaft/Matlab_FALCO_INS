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
    % КАЧЕСТВО ОЦЕНКИ BIAS — строго в момент ПЕРЕД терминальным коастом.
    % Сравнивается оценка фильтра с ПОЛНЫМ истинным смещением в тот же
    % момент: turn-on + текущая in-run составляющая (+ VRC, который фильтр
    % от bias не отличает). Сравнение с одним лишь turn-on некорректно:
    % у гироскопа Pulse-40 in-run bias РАВЕН калиброванному turn-on,
    % отчего метрика "сколько выучено" давала -536%.
    ba_true_pre  = nan(N_mc,3);   ba_est_pre  = nan(N_mc,3);
    bg_true_pre  = nan(N_mc,3);   bg_est_pre  = nan(N_mc,3);
    ba_resid_pre = nan(N_mc,1);   bg_resid_pre = nan(N_mc,1);
    % МЕТРИКИ СОСТОЯТЕЛЬНОСТИ ФИЛЬТРА (NEES / NIS).
    % Для СОГЛАСОВАННОГО фильтра среднее обеих метрик равно размерности:
    %   NEES6 (позиция+скорость) -> 6
    %   NIS   (измерение GNSS)   -> 6
    % Систематическое превышение означает ПЕРЕУВЕРЕННОСТЬ: реальные ошибки
    % больше, чем фильтр ожидает по своей ковариации. Это ожидаемо, поскольку
    % в модели истины есть источники, которых нет в состоянии фильтра
    % (перекос осей, g-sens, B2, B4).
    % ХРАНИМ ПОЛНЫЕ ВРЕМЕННЫЕ РЯДЫ, а не усреднённые по времени скаляры.
    % Причина: отсчёты NEES/NIS ВНУТРИ одной реализации КОРРЕЛИРОВАНЫ
    % (состояние фильтра меняется плавно), поэтому их среднее по времени
    % НЕЛЬЗЯ трактовать как выборку независимых chi-квадрат величин.
    % Корректная статистика — ANEES/ANIS: усреднение по АНСАМБЛЮ при
    % ФИКСИРОВАННОМ времени, где реализации независимы по построению.
    NEES_mat = [];              % (N_mc x n_diag)  заполняется по первой реализации
    NIS_mat  = [];              % (N_mc x n_gnss)
    t_diag   = [];              % общая диагностическая сетка
    t_gnss   = [];              % общая сетка обновлений GNSS

    % Усреднённые по времени величины — ТОЛЬКО как компактная сводка.
    % Строгим 95% тестом они НЕ являются (см. выше про корреляцию).
    nees_gnss  = nan(N_mc,1);
    nees_coast = nan(N_mc,1);
    nis_tavg   = nan(N_mc,1);

    % Состояние непосредственно перед терминальным коастом
    pre_horiz  = nan(N_mc,1);
    pre_dv     = nan(N_mc,1);
    pre_dpsi   = nan(N_mc,1);
    align_used = zeros(N_mc,3);
    % Применённые начальные ошибки (диагностика)
    init_pos_used = zeros(N_mc,3);
    init_vel_used = zeros(N_mc,3);
    n_fail     = 0;

    % Диагностика фактически разыгранных аномалий поля
    dg_info = struct('defl_E_arcsec',cell(N_mc,1), 'defl_N_arcsec',[], ...
                     'defl_total_arcsec',[], 'anom_vert_mgal',[], ...
                     'dg_horiz_mps2',[], 'dg_norm_mps2',[]);
    dg_vec  = zeros(N_mc,3);

    % Начало терминального коаста: второе окно аутэйджа. Метрики "перед
    % коастом" берутся там, где фильтр в последний раз имел GNSS —
    % это разделяет вклад активного участка и бесспутникового.
    if size(cfg.outage,1) >= 2
        t_coast_start = cfg.outage(2,1);
    else
        t_coast_start = truth.t(end);
    end

    for i = 1:N_mc
        seed_i = base_seed + i;

        % --- Розыгрыш ошибок по профилю ИСТИНЫ ---
        rng_imu  = RandStream('mt19937ar','Seed', seed_i);
        e_imu    = imu_draw_errors(prof_truth, rng_imu);
        imu_meas = add_imu_errors(imu_ideal, prof_truth, e_imu, rng_imu);

        % --- Ошибка выставки (своя на каждой реализации) ---
        cfg_i = cfg;
        cfg_i.align_err = cfg.align_sigma(:) .* randn(rng_imu, 3, 1);

        % НАЧАЛЬНАЯ ОШИБКА ПОЗИЦИИ И СКОРОСТИ (предстартовое усреднение GNSS).
        % Разыгрывается из ТОГО ЖЕ потока rng_imu сразу после выставки, чтобы
        % последовательность оставалась воспроизводимой при фиксированном
        % base_seed. Каждая реализация получает свою начальную ошибку.
        % При отсутствии полей в cfg ошибка нулевая (идеальный старт).
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
            ba_resid_pre(i) = NaN;
            bg_resid_pre(i) = NaN;
        else
            err_horiz(i) = res.err_horiz_final;
            err_3d(i)    = res.err_3d_final;
            err_enu(i,:) = res.dr_enu(end,:);
            % Строго ПОСЛЕДНЯЯ диагностическая точка ДО начала аутэйджа
            k_pre = find(res.t < t_coast_start, 1, 'last');
            if isempty(k_pre), k_pre = 1; end

            % --- Метрики состоятельности: СОХРАНЯЕМ РЯДЫ ---
            if isfield(res,'nees6')
                if isempty(NEES_mat)
                    t_diag   = res.t(:)';
                    NEES_mat = nan(N_mc, numel(t_diag));
                end
                % Все реализации идут по одной траектории с одним diag_decim,
                % поэтому сетка совпадает. Проверка на случай расхождения.
                if numel(res.nees6) == size(NEES_mat,2)
                    NEES_mat(i,:) = res.nees6(:)';
                end

                % Компактная сводка по участкам (НЕ строгий тест)
                i_coast = res.t >= t_coast_start;
                i_gnss  = ~i_coast & (res.t > cfg.outage(1,2));
                if any(i_gnss),  nees_gnss(i)  = mean(res.nees6(i_gnss),  'omitnan'); end
                if any(i_coast), nees_coast(i) = mean(res.nees6(i_coast), 'omitnan'); end
            end
            if isfield(res,'gnss_nis') && ~isempty(res.gnss_nis)
                if isempty(NIS_mat)
                    t_gnss  = res.gnss_t(:)';
                    NIS_mat = nan(N_mc, numel(t_gnss));
                end
                if numel(res.gnss_nis) == size(NIS_mat,2)
                    NIS_mat(i,:) = res.gnss_nis(:)';
                end
                nis_tavg(i) = mean(res.gnss_nis, 'omitnan');
            end

            pre_horiz(i) = hypot(res.dr_ned(k_pre,1), res.dr_ned(k_pre,2));
            pre_dv(i)    = norm(res.dv(k_pre,:));
            pre_dpsi(i)  = norm(res.dpsi(k_pre,:)) * 1e3;   % мрад

            % Истинное ПОЛНОЕ смещение в тот же момент: берётся из истории,
            % сохранённой add_imu_errors, а не пересчитывается заново.
            [~, k_truth] = min(abs(imu_meas.t - res.t(k_pre)));
            bg_true_pre(i,:) = imu_meas.bg_total(k_truth,:);
            ba_true_pre(i,:) = imu_meas.ba_total(k_truth,:);
            bg_est_pre(i,:)  = res.bg_est(k_pre,:);
            ba_est_pre(i,:)  = res.ba_est(k_pre,:);

            bg_resid_pre(i) = norm(bg_true_pre(i,:) - bg_est_pre(i,:));
            ba_resid_pre(i) = norm(ba_true_pre(i,:) - ba_est_pre(i,:));
        end

        align_used(i,:)     = cfg_i.align_err';
        init_pos_used(i,:)  = cfg_i.init_pos_err';
        init_vel_used(i,:)  = cfg_i.init_vel_err';
    end

    % --- Статистика ---
    ok = ~isnan(err_horiz);
    mc.err_horiz  = err_horiz;    mc.err_3d  = err_3d;    mc.err_enu = err_enu;
    mc.align_used     = align_used;
    mc.init_pos_used  = init_pos_used;
    mc.init_vel_used  = init_vel_used;
    % Качество оценки bias перед коастом (векторы и остатки)
    mc.ba_true_pre  = ba_true_pre;    mc.ba_est_pre  = ba_est_pre;
    mc.bg_true_pre  = bg_true_pre;    mc.bg_est_pre  = bg_est_pre;
    mc.ba_resid_pre = ba_resid_pre;   mc.bg_resid_pre = bg_resid_pre;
    mc.pre_horiz  = pre_horiz;    mc.pre_dv  = pre_dv;    mc.pre_dpsi = pre_dpsi;
    mc.t_coast_start = t_coast_start;
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

    % Медианы метрик перед коастом
    % =====================================================================
    % СТАТИСТИКА СОСТОЯТЕЛЬНОСТИ: ANEES / ANIS
    % =====================================================================
    % ANEES(k) = среднее NEES по ансамблю в МОМЕНТ t(k).
    % Реализации независимы, поэтому при согласованном фильтре
    %       sum_i NEES_i(k) ~ chi2(N*dim)
    % и доверительный коридор для ANEES(k):
    %       [ chi2(0.025, N*dim)/N , chi2(0.975, N*dim)/N ]
    % Квантили — аппроксимация Уилсона-Хилферти (без Statistics Toolbox).
    %
    % Так усреднение идёт ПОПЕРЁК независимых реализаций, а не вдоль
    % коррелированного времени, поэтому тест корректен.
    dim = 6;
    mc.consist_dim = dim;

    if ~isempty(NEES_mat)
        n_ok_k = sum(~isnan(NEES_mat), 1);              % реализаций в каждый момент
        mc.anees   = mean(NEES_mat, 1, 'omitnan');      % ANEES(t)
        mc.anees_t = t_diag;
        mc.anees_n = n_ok_k;
        % Коридор строится по фактическому числу реализаций в момент k
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
        % Типичное число реализаций в момент времени (для подписи графика)
        mc.consist_n_typ = median(n_ok_k(n_ok_k > 0));
    end

    if ~isempty(NIS_mat)
        n_ok_g = sum(~isnan(NIS_mat), 1);
        mc.anis   = mean(NIS_mat, 1, 'omitnan');        % ANIS(t)
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

    % Компактные средние по времени. ВНИМАНИЕ: это НЕ строгий 95% тест,
    % поскольку отсчёты внутри реализации коррелированы. Служат только
    % грубым индикатором уровня.
    mc.nees_gnss  = nees_gnss;
    mc.nees_coast = nees_coast;
    mc.nis_tavg   = nis_tavg;
    mc.nees_gnss_mean  = mean(nees_gnss(~isnan(nees_gnss)));
    mc.nees_coast_mean = mean(nees_coast(~isnan(nees_coast)));
    mc.nis_tavg_mean   = mean(nis_tavg(~isnan(nis_tavg)));

    mc.pre_horiz_med = median(pre_horiz(~isnan(pre_horiz)));
    mc.pre_dv_med    = median(pre_dv(~isnan(pre_dv)));
    mc.pre_dpsi_med  = median(pre_dpsi(~isnan(pre_dpsi)));

    % Медианы истинного смещения и остатка оценки перед коастом.
    % ПРОЦЕНТ "сколько выучено" СОЗНАТЕЛЬНО НЕ считается: при малом
    % истинном смещении отношение неустойчиво и даёт бессмысленные
    % значения. Смотреть следует на пару (истинное, остаток).
    mc.bg_true_pre_med  = median(vecnorm(bg_true_pre(~isnan(bg_resid_pre),:),2,2));
    mc.ba_true_pre_med  = median(vecnorm(ba_true_pre(~isnan(ba_resid_pre),:),2,2));
    mc.bg_resid_pre_med = median(bg_resid_pre(~isnan(bg_resid_pre)));
    mc.ba_resid_pre_med = median(ba_resid_pre(~isnan(ba_resid_pre)));
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
    % Для conf = 0.95 это ±1.959964; для прочих значений используется
    % обратная функция ошибок (входит в базовый MATLAB).
    z = sqrt(2) * erfinv(1 - alpha);

    lo = k * (1 - 2/(9*k) - z*sqrt(2/(9*k)))^3;
    hi = k * (1 - 2/(9*k) + z*sqrt(2/(9*k)))^3;
end
