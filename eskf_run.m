function res = eskf_run(imu_meas, truth, c, prof, cfg)
%ESKF_RUN  Главный цикл ESKF: механизация + predict + GNSS update + feedback.
%
%   Полный контур интегрированной навигации ИНС/GNSS в ECEF с моделью
%   пропаданий GNSS (аутэйдж на разгоне и на терминальном коасте).
%
%   Цикл на каждом такте навигатора:
%     1. Получить приращения IMU (dtheta, dvel) за такт
%     2. Механизация: скорректировать измерения, проинтегрировать (r,v,q)
%     3. Predict: построить F -> PHI, Qd -> распространить ковариацию P
%     4. Если GNSS доступен: построить H,R -> невязка -> update -> feedback+reset
%     5. Сохранить диагностику
%
%   Вход:
%     imu_meas - ЗАШУМЛЁННЫЕ показания IMU (.t, .fb, .wib_b, .Ceb)
%                ОБЯЗАТЕЛЬНО .w_eb_b — истинная угловая скорость тела
%                относительно ECEF (служебная кинематика истины, НЕ измерение).
%                Нужна для синтеза истинной скорости GNSS-антенны. Проброс
%                выполняет add_imu_errors из attitude_program.
%     truth    - истинная траектория (.t, .R, .V, .dt, .r0, .Cen)
%                опционально .C_e_ned0 (NED -> ECEF). Если поля нет (старые
%                .mat), матрица выводится из .Cen перестановкой столбцов.
%     c        - константы
%     prof     - профиль IMU (для Qc и начальной ковариации)
%     cfg      - конфигурация прогона (см. main_step5.m)
%                опционально .init_pos_err / .init_vel_err (3x1) — начальная
%                ошибка позиции/скорости навигатора. При отсутствии полей
%                инициализация ИДЕАЛЬНАЯ (нулевые ошибки), что сохраняет
%                поведение детерминированных тестов.
%   Выход:
%     res      - структура результатов и диагностики

    % =====================================================================
    % Шаг 5.11.1: ПОДГОТОВКА - приращения IMU на такте навигатора
    % =====================================================================
    [a_idx, dtheta, dvel, dt, step] = make_increments(imu_meas, truth, cfg.fnav);
    M = numel(a_idx);

    if ~isfield(cfg,'lever')      || isempty(cfg.lever),      cfg.lever = zeros(3,1); end
    if ~isfield(cfg,'diag_decim') || isempty(cfg.diag_decim), cfg.diag_decim = 25;    end
    rng_gnss = RandStream('mt19937ar','Seed', cfg.seed);

    % Проверка наличия служебной истинной кинематики. Fallback на зашумлённый
    % imu_meas.wib_b СОЗНАТЕЛЬНО не делается: истина не должна строиться из
    % показаний датчика. Отсутствие поля означает устаревший pipeline.
    assert(isfield(imu_meas, 'w_eb_b'), 'eskf_run:missingWebB', ...
        ['imu_meas.w_eb_b отсутствует. Это истинная угловая скорость тела\n' ...
         'относительно ECEF, нужная для синтеза скорости GNSS-антенны.\n' ...
         'Поле формируется в attitude_program и пробрасывается в\n' ...
         'add_imu_errors. Перегенерируй измерения актуальной версией.']);

    % =====================================================================
    % Шаг 5.11.2: ИНИЦИАЛИЗАЦИЯ номинального состояния (с ошибкой выставки)
    % =====================================================================
    C_true0 = imu_meas.Ceb(:,:,1);
    % C_true = (I + [dpsi x]) * C_hat  =>  C_hat = (I - [dpsi x]) * C_true
    C_hat0  = (eye(3) - skew(cfg.align_err)) * C_true0;

    % =====================================================================
    % НАЧАЛЬНАЯ ОШИБКА ПОЗИЦИИ И СКОРОСТИ
    % =====================================================================
    % Физически это результат предстартового усреднения решения GNSS на
    % неподвижном изделии: позиция и скорость известны не точно, а с
    % остаточной погрешностью. Прежде навигатор стартовал ТОЧНО из истины,
    % что занижало ошибку на начальном участке.
    %
    % ОБРАТНАЯ СОВМЕСТИМОСТЬ: при отсутствии полей ошибки нулевые, то есть
    % инициализация остаётся идеальной. Это сохраняет поведение
    % детерминированных тестов (numerical floor, gravity anomaly,
    % gnss time offset), которые этих полей не задают.
    if isfield(cfg, 'init_pos_err')
        init_pos_err = cfg.init_pos_err(:);
    else
        init_pos_err = zeros(3,1);
    end
    if isfield(cfg, 'init_vel_err')
        init_vel_err = cfg.init_vel_err(:);
    else
        init_vel_err = zeros(3,1);
    end

    nav.r  = truth.R(1,:)' - init_pos_err;  % Знак "-" выбран специально для соответствия принятой конвенции ESKF. 
    nav.v  = truth.V(1,:)' - init_vel_err;
    nav.q  = C2q(C_hat0);
    nav.C  = q2C(nav.q);
    nav.bg = zeros(3,1);
    nav.ba = zeros(3,1);
    nav.sg = zeros(3,1);
    nav.sa = zeros(3,1);
    nav.dth_prev = zeros(3,1);
    nav.dv_prev  = zeros(3,1);

    % Паспортная g-чувствительность гироскопа для компенсации в механизации.
    % Берётся из imu_meas (её кладёт туда add_imu_errors) либо из профиля.
    % Без неё на траектории с ненулевой удельной силой на коасте копится
    % ~2 мрад ошибки ориентации за 60 с и ~8 м систематического бокового сноса.
    if isfield(imu_meas, 'Gg_nom')
        nav.Gg_cal = imu_meas.Gg_nom;
    elseif isfield(prof.gyro, 'g_sensitivity')
        nav.Gg_cal = prof.gyro.g_sensitivity * eye(3);
    else
        nav.Gg_cal = zeros(3);
    end

    % Ошибка гравитационной модели навигатора (аномалия поля).
    % Передаётся снаружи (разыграна на реализацию МК) либо ноль.
    if isfield(cfg, 'dg_model')
        nav.dg_model = cfg.dg_model;
    else
        nav.dg_model = zeros(3,1);
    end

    % =====================================================================
    % Шаг 5.11.3: ИНИЦИАЛИЗАЦИЯ фильтра (P0)
    % =====================================================================
    kf.x = zeros(21,1);
    kf.P = zeros(21,21);
    kf.P(1:3,1:3)     = diag(cfg.P0_pos(:).^2);
    kf.P(4:6,4:6)     = diag(cfg.P0_vel(:).^2);
    kf.P(7:9,7:9)     = diag(cfg.P0_att(:).^2);
    % =====================================================================
    % НАЧАЛЬНАЯ КОВАРИАЦИЯ СМЕЩЕНИЙ — СОГЛАСОВАНА С ИСТИНОЙ
    % =====================================================================
    % ИСПРАВЛЕНО. Прежде здесь стояло turnon_bias_sigma^2 БЕЗ учёта
    % предпусковой калибровки, тогда как истина формирует смещение так:
    %
    %   imu_draw_errors : turnon = cal_factor * turnon_bias_sigma * randn
    %   add_imu_errors  : b_inrun(0) = inrun_bias_sigma * randn
    %   полное смещение : bg_total = turnon + b_inrun + vrc
    %
    % Фактическое стартовое СКО складывается из ДВУХ независимых частей:
    %
    %   sigma_b0 = sqrt( (cal_factor * turnon_bias_sigma)^2
    %                    + inrun_bias_sigma^2 )
    %
    % Для гироскопа Pulse-40 это sqrt(0.495^2 + 0.5^2) = 0.704 °/ч на ось,
    % тогда как в фильтр подставлялось 150 °/ч — завышение в 213 раз по СКО
    % и в 45 000 раз по дисперсии. Настолько раздутая P0 давала фильтру
    % «право» списывать на bias посторонние эффекты (остаток g-sens,
    % перекос осей, ошибки масштаба), из-за чего оценка уходила на 3-4 °/ч
    % по осям X и Z, где смещение почти не наблюдаемо.
    %
    % Та же ошибка была у акселерометра, но там завышение всего вчетверо,
    % поэтому проявлялась она слабее.
    %
    % VRC в P0 НЕ включается: это не начальная неопределённость, а
    % детерминированный отклик на вибрацию, действующий только при
    % работающем двигателе.
    cal_g = get_cal(prof, 'cal_factor_gyro');
    cal_a = get_cal(prof, 'cal_factor_accel');

    sig_bg0 = sqrt( (cal_g * prof.gyro.turnon_bias_sigma)^2 ...
                    + prof.gyro.inrun_bias_sigma^2 );
    sig_ba0 = sqrt( (cal_a * prof.accel.turnon_bias_sigma)^2 ...
                    + prof.accel.inrun_bias_sigma^2 );

    kf.P(10:12,10:12) = sig_bg0^2 * eye(3);
    kf.P(13:15,13:15) = sig_ba0^2 * eye(3);
    kf.P(16:18,16:18) = prof.gyro.scale_factor_sigma^2  * eye(3);
    kf.P(19:21,19:21) = prof.accel.scale_factor_sigma^2 * eye(3);

    % =====================================================================
    % Шаг 5.11.4: ПРЕАЛЛОКАЦИЯ диагностики
    % =====================================================================
    n_diag = ceil(M/cfg.diag_decim) + 1;
    d_t    = zeros(n_diag,1);
    d_dr   = zeros(n_diag,3);   d_dv    = zeros(n_diag,3);   d_dpsi  = zeros(n_diag,3); % Ошибки поз, скорость, ориентация. 
    d_sr   = zeros(n_diag,3);   d_sv    = zeros(n_diag,3);   d_sp    = zeros(n_diag,3);
    d_ba   = zeros(n_diag,3);   d_bg    = zeros(n_diag,3);
    d_sba  = zeros(n_diag,3);   d_sbg   = zeros(n_diag,3);
    % Оценки МАСШТАБНЫХ КОЭФФИЦИЕНТОВ и их СКО.
    % В текущей модели scale — случайная константа (F = 0, Qc = 0), поэтому
    % P по этим состояниям может только УБЫВАТЬ: рост неопределённости не
    % предусмотрен. Диагностика показывает, насколько фильтр фактически
    % использует эту возможность при имеющемся профиле возбуждения.
    d_sg   = zeros(n_diag,3);   d_sa    = zeros(n_diag,3);
    d_ssg  = zeros(n_diag,3);   d_ssa   = zeros(n_diag,3);
    % ИСТИННОЕ ПОЛНОЕ смещение гироскопа (turn-on + in-run + VRC) — из
    % служебной истории add_imu_errors. Нужно потому, что фильтр оценивает
    % ПОЛНОЕ смещение, и сравнивать его только с turn-on некорректно:
    % у Pulse-40 in-run bias равен калиброванному turn-on.
    d_bgtr = zeros(n_diag,3);
    % КОРРЕЛЯЦИЯ δbg_i <-> δsg_i по каждой оси.
    % Модель измерения гироскопа: δω = δbg + diag(ω)·δsg. При слабо
    % меняющейся ω состояния почти НЕРАЗДЕЛИМЫ: |rho| -> 1 означает, что
    % фильтр не может развести их, и ошибка одного гасится другим.
    % Это прямая мера разделимости, и она уже содержится в самой P.
    d_rbgsg= zeros(n_diag,3);
    % Абсолютная навигационная траектория (не ошибка, а само решение)
    d_rnav = zeros(n_diag,3);   d_vnav  = zeros(n_diag,3);
    % NEES по позиции и скорости (усечённый, 6 состояний).
    % Полный NEES по 21 состоянию потребовал бы истинных ошибок bias и scale,
    % которые здесь недоступны; усечённый вариант использует только dr и dv,
    % а для них истина известна точно.
    d_nees6 = zeros(n_diag,1);
    % РАЗДЕЛЬНЫЕ NEES по позиции и скорости (dim = 3 каждый).
    % Полный NEES6 смешивает два канала; раздельные показывают, в каком
    % именно из них фильтр переуверен — это разные физические причины.
    d_neesr = zeros(n_diag,1);
    d_neesv = zeros(n_diag,1);
    id = 0;

    % ОТДЕЛЬНЫЙ лог невязок GNSS (пишется В МОМЕНТ обновления, не по прореживанию!)
    % Иначе такты диагностики и такты GNSS не совпадают -> лог остаётся пустым.
    n_gnss_max = ceil((truth.t(end) - truth.t(1)) * cfg.f_gnss) + 2;
    g_t     = zeros(n_gnss_max,1);
    g_innov = zeros(n_gnss_max,6);
    % Диагностика временного сдвига (B4)
    g_nis      = zeros(n_gnss_max,1);   % NIS = z'·S^-1·z (состоятельность)
    % Отношение diag(H·P·H') / diag(R) по каждому обновлению.
    % Проверяет, насколько шум измерений доминирует в инновационной
    % ковариации S = H·P·H' + R. При отношении << 1 величина S почти
    % целиком определяется R, и NIS теряет чувствительность к
    % заниженности P.
    g_rho_pos  = zeros(n_gnss_max,3);   % компоненты 1:3 (позиция)
    g_rho_vel  = zeros(n_gnss_max,3);   % компоненты 4:6 (скорость)
    g_t_meas   = zeros(n_gnss_max,1);   % момент, к которому относится измерение
    g_age      = zeros(n_gnss_max,1);   % возраст измерения t_nav - t_meas
    g_mismatch = zeros(n_gnss_max,3);   % r_ant_true(t_meas) - r_ant_true(t_nav)
    ig = 0;

    t_next_gnss = truth.t(1);
    dt_gnss     = 1/cfg.f_gnss;

    % =====================================================================
    % Шаг 5.11.4b: ВРЕМЕННОЙ СДВИГ GNSS (B4)
    % =====================================================================
    % Моделируется НЕКОМПЕНСИРОВАННАЯ задержка измерения:
    %     t_meas = t_nav - cfg.gnss_time_offset
    % Измерение синтезируется по истине в момент t_meas, а коррекция фильтра
    % выполняется в текущий момент t_nav — предсказание, H и update остаются
    % привязанными к t_nav. Именно в этом и состоит эффект рассинхронизации.
    %
    % Сдвиг реализуется ЦЕЛОЧИСЛЕННЫМ смещением индекса сетки истины.
    % Интерполяция сознательно не применяется: при типовых сдвигах
    % 0.5..10 мс и сетке 2000 Гц (шаг 0.5 мс) все значения кратны шагу,
    % а погрешность интерполяции при 0.5 мс могла бы превысить сам
    % измеряемый эффект. Некратный сдвиг прерывается assert-ом.
    if isfield(cfg,'gnss_time_offset_enable') && cfg.gnss_time_offset_enable
        gnss_dt_off = cfg.gnss_time_offset;
    else
        gnss_dt_off = 0.0;
    end

    % Знак: измерение может быть только УСТАРЕВШИМ, не из будущего
    assert(gnss_dt_off >= 0, 'eskf_run:negativeGnssOffset', ...
        ['cfg.gnss_time_offset = %.6g с < 0. Отрицательный сдвиг означал бы\n' ...
         'измерение из будущего, что физически невозможно.'], gnss_dt_off);

    n_shift_raw = gnss_dt_off / truth.dt;
    n_shift     = round(n_shift_raw);

    % Кратность шагу сетки истины
    assert(abs(n_shift_raw - n_shift) <= 1e-9 * max(1, n_shift_raw), ...
        'eskf_run:gnssOffsetNotCommensurate', ...
        ['cfg.gnss_time_offset = %.6g с НЕ КРАТЕН шагу сетки истины %.6g с.\n' ...
         'Отношение = %.9f, ближайшее целое = %d.\n' ...
         'Дробный сдвиг потребовал бы интерполяции (для Ceb — slerp),\n' ...
         'что на данном этапе не реализовано.'], ...
        gnss_dt_off, truth.dt, n_shift_raw, n_shift);

    % =====================================================================
    % Шаг 5.11.5: ГЛАВНЫЙ ЦИКЛ
    % =====================================================================
    for k = 1:M
        a  = a_idx(k);
        b  = a + step;
        tk = truth.t(b);

        % --- (1) МЕХАНИЗАЦИЯ: один такт ---
        [nav, aux] = ins_mech_step(nav, dtheta(k,:)', dvel(k,:)', dt, c);

        % --- (2) PREDICT: F -> PHI -> Qd -> P ---
        F   = eskf_F_matrix(nav.r, nav.v, nav.C, aux.f_b, aux.w_b, c, cfg.corrtime);
        PHI = eskf_discretize(F, dt, 2);
        Qd  = eskf_Q_matrix(nav.C, PHI, dt, prof, cfg.corrtime);
        kf.P = PHI * kf.P * PHI' + Qd;
        kf.P = (kf.P + kf.P')/2;

        % --- (3) GNSS UPDATE (если пришло время И нет аутэйджа) ---
        gnss_used = false;
        if tk >= t_next_gnss
            % Индекс сетки истины, к которому относится измерение.
            % При n_shift = 0 совпадает с b (идеальная синхронизация).
            b_meas = b - n_shift;

            % Измерение доступно только если сдвинутый момент лежит внутри
            % траектории. В первые миллисекунды полёта это не так —
            % обновление ПРОПУСКАЕТСЯ (индекс не зажимается, иначе возник бы
            % фиктивный сдвиг, отличный от заданного).
            if gnss_available(tk, cfg.outage) && b_meas >= 1
                % ---- ИСТИННЫЕ позиция и скорость АНТЕННЫ ----
                % Измерение и предсказание ОБЯЗАНЫ относиться к одной
                % физической точке. Раньше измерение строилось в точке IMU,
                % а предсказание — в точке антенны, что давало систематику
                % порядка |lever| (наблюдалось 0.155 м при |lever| = 0.141 м).
                % ВСЯ ЧЕТВЁРКА берётся СОГЛАСОВАННО в момент t_meas (индекс
                % b_meas). Сдвигать только R/V, оставив Ceb в t_nav, нельзя:
                % это внесло бы искусственную ошибку ~|w x l|*dt, которой в
                % действительности нет.
                C_true_b   = imu_meas.Ceb(:,:,b_meas);     % истинная ориентация
                w_eb_true  = imu_meas.w_eb_b(b_meas,:)';   % истинная w отн. ECEF
                r_ant_true = truth.R(b_meas,:)' + C_true_b * cfg.lever;
                v_ant_true = truth.V(b_meas,:)' ...
                             + C_true_b * cross(w_eb_true, cfg.lever);

                % Синтез измерения GNSS: истина АНТЕННЫ + шум
                r_gnss = r_ant_true + cfg.sigma_pos .* randn(rng_gnss,3,1);
                v_gnss = v_ant_true + cfg.sigma_vel .* randn(rng_gnss,3,1);

                % ---- Предсказание: позиция/скорость АНТЕННЫ по оценке ----
                % Земной член входит со знаком МИНУС:
                %   v_ant = v + C(w_ib x l) - [w_ie x](C l),
                % так как w_eb^b = w_ib^b - C'*w_ie и
                %   C*((C'w_ie) x l) = (C C' w_ie) x (C l) = w_ie x (C l).
                % Ранее стоял плюс; расхождение 2*|w_ie x C*l| ~ 2e-5 м/с.
                r_ant_hat = nav.r + nav.C * cfg.lever;
                v_ant_hat = nav.v + nav.C * cross(aux.w_b, cfg.lever) ...
                                  - skew(c.WIE) * (nav.C * cfg.lever);

                % НЕВЯЗКА: ИЗМЕРЕНИЕ минус оценка (конвенция Sola)
                z = [r_gnss - r_ant_hat; v_gnss - v_ant_hat];

                [H, R] = eskf_H_gnss(cfg.sigma_pos, cfg.sigma_vel, true, ...
                                     nav.C, aux.w_b, cfg.lever, c);
                kf = eskf_update(kf, z, H, R);

                % NIS — нормированный квадрат невязки. Для состоятельного
                % фильтра E[NIS] = dim(z) = 6. Значения СИСТЕМАТИЧЕСКИ выше
                % означают, что фильтр ПЕРЕУВЕРЕН: реальные невязки больше,
                % чем он ожидает по своей ковариации.
                nis_k = z' * (kf.S \ z);

                % Вклад P против вклада R в ИННОВАЦИОННОЙ ковариации.
                % ИСПРАВЛЕНО: раньше бралось diag(H*kf.P*H') ПОСЛЕ
                % eskf_update, то есть по АПОСТЕРИОРНОЙ P, тогда как NIS
                % использует S = H*P_prior*H' + R. Матрицы были разные, и
                % отношение не соответствовало тому, что видит NIS.
                % Теперь H*P_prior*H' извлекается из той же S: S - R.
                HPH_diag = diag(kf.S) - diag(R);
                R_diag   = diag(R);
                rho_k    = (HPH_diag ./ max(R_diag, eps))';

                [kf, nav] = eskf_feedback(kf, nav);

                gnss_used = true;

                % Лог невязки СРАЗУ (независимо от прореживания диагностики)
                ig = ig + 1;
                g_t(ig)       = tk;
                g_innov(ig,:) = z';
                g_nis(ig)     = nis_k;
                % Каналы НЕ смешиваются: масштабы R по позиции (м) и по
                % скорости (м/с) различаются на два порядка, поэтому
                % усреднение шести компонент в одно число скрывало бы
                % картину.
                g_rho_pos(ig,:) = rho_k(1:3);
                g_rho_vel(ig,:) = rho_k(4:6);

                % --- Диагностика временного сдвига (B4) ---
                g_t_meas(ig) = truth.t(b_meas);
                g_age(ig)    = tk - truth.t(b_meas);
                % Рассогласование истины: где антенна БЫЛА в t_meas против
                % того, где она находится в t_nav. Это чистый геометрический
                % эффект задержки, не зависящий от работы фильтра.
                C_true_now      = imu_meas.Ceb(:,:,b);
                r_ant_true_now  = truth.R(b,:)' + C_true_now * cfg.lever;
                g_mismatch(ig,:) = (r_ant_true - r_ant_true_now)';
            end
            t_next_gnss = t_next_gnss + dt_gnss;
        end

        % --- (4) ДИАГНОСТИКА (с прореживанием) ---
        if mod(k-1, cfg.diag_decim) == 0 || k == M
            id = id + 1;
            d_t(id)     = tk;
            d_dr(id,:)  = (nav.r - truth.R(b,:)')';
            d_dv(id,:)  = (nav.v - truth.V(b,:)')';
            dC          = imu_meas.Ceb(:,:,b) * nav.C' - eye(3);
            d_dpsi(id,:)= [dC(3,2), dC(1,3), dC(2,1)];
            sP          = sqrt(diag(kf.P));
            d_sr(id,:)  = sP(1:3)';    d_sv(id,:)  = sP(4:6)';    d_sp(id,:) = sP(7:9)';
            d_bg(id,:)  = nav.bg';     d_ba(id,:)  = nav.ba';
            d_sbg(id,:) = sP(10:12)';  d_sba(id,:) = sP(13:15)';
            % Накопленная оценка scale после feedback/injection
            d_sg(id,:)  = nav.sg';     d_sa(id,:)  = nav.sa';
            d_ssg(id,:) = sP(16:18)';  d_ssa(id,:) = sP(19:21)';

            % Истинное полное смещение гироскопа в этот момент
            if isfield(imu_meas,'bg_total')
                d_bgtr(id,:) = imu_meas.bg_total(b,:);
            end

            % Корреляция δbg_i <-> δsg_i: индексы состояний 9+i и 15+i
            for iax = 1:3
                ib  = 9  + iax;
                isc = 15 + iax;
                den = sqrt(kf.P(ib,ib) * kf.P(isc,isc));
                if den > 0
                    d_rbgsg(id,iax) = kf.P(ib,isc) / den;
                end
            end
            % Абсолютное навигационное решение в ECEF
            d_rnav(id,:) = nav.r';     d_vnav(id,:) = nav.v';

            % NEES по блоку позиция/скорость: dx' · P^-1 · dx, где dx —
            % ИСТИННАЯ ошибка состояния (истина минус оценка, конвенция Sola).
            % Для состоятельного фильтра E[NEES6] = 6.
            % P берётся ДО следующего предсказания, согласованно с dx.
            dr_true = truth.R(b,:)' - nav.r;
            dv_true = truth.V(b,:)' - nav.v;
            dx6   = [dr_true; dv_true];
            P6    = kf.P(1:6,1:6);
            P6    = (P6 + P6')/2;                 % симметризация
            d_nees6(id) = dx6' * (P6 \ dx6);

            % Раздельно по каналам: берутся ДИАГОНАЛЬНЫЕ блоки Prr и Pvv.
            % Кросс-корреляция позиция-скорость при этом отбрасывается,
            % поэтому сумма NEES_r + NEES_v в общем случае НЕ равна NEES6.
            Prr = kf.P(1:3,1:3);   Prr = (Prr + Prr')/2;
            Pvv = kf.P(4:6,4:6);   Pvv = (Pvv + Pvv')/2;
            d_neesr(id) = dr_true' * (Prr \ dr_true);
            d_neesv(id) = dv_true' * (Pvv \ dv_true);
        end
    end

    % =====================================================================
    % Шаг 5.11.6: УПАКОВКА РЕЗУЛЬТАТА
    % =====================================================================
    idx = 1:id;
    res.t      = d_t(idx);
    res.dr     = d_dr(idx,:);     res.dv    = d_dv(idx,:);    res.dpsi  = d_dpsi(idx,:);
    res.sig_r  = d_sr(idx,:);     res.sig_v = d_sv(idx,:);    res.sig_p = d_sp(idx,:);
    res.bg_est = d_bg(idx,:);     res.ba_est= d_ba(idx,:);
    res.sig_bg = d_sbg(idx,:);    res.sig_ba= d_sba(idx,:);
    % Масштабные коэффициенты (безразмерные; в ppm умножать на 1e6)
    res.sg_est = d_sg(idx,:);     res.sa_est = d_sa(idx,:);
    res.sig_sg = d_ssg(idx,:);    res.sig_sa = d_ssa(idx,:);
    % Диагностика смещения гироскопа
    res.bg_true   = d_bgtr(idx,:);      % полное истинное смещение [рад/с]
    res.rho_bg_sg = d_rbgsg(idx,:);     % корреляция δbg_i <-> δsg_i

    % Абсолютная навигационная траектория (ECEF). Логируется всегда:
    % нужна для построения truth vs nav и для внешнего NED-интерфейса.
    res.r_nav  = d_rnav(idx,:);
    res.v_nav  = d_vnav(idx,:);

    % --- Метрики состоятельности фильтра ---
    % NEES6: E = 6 для состоятельного фильтра. Больше -> переуверен.
    % NIS:   E = 6 (размерность измерения GNSS). Больше -> переуверен.
    res.nees6      = d_nees6(idx);
    res.nees6_dim  = 6;
    res.nees_r     = d_neesr(idx);      % dim = 3
    res.nees_v     = d_neesv(idx);      % dim = 3
    res.nees_rv_dim= 3;
    res.nis_dim    = 6;

    % Невязки GNSS - в отдельных полях со своей временной шкалой
    gidx           = 1:ig;
    res.gnss_t     = g_t(gidx);
    res.gnss_innov = g_innov(gidx,:);
    res.n_gnss     = ig;

    % --- Диагностика временного сдвига GNSS (B4) ---
    res.gnss_nis       = g_nis(gidx);         % NIS по каждому обновлению
    res.gnss_rho_pos   = g_rho_pos(gidx,:);   % diag(HPH')/diag(R), позиция
    res.gnss_rho_vel   = g_rho_vel(gidx,:);   % то же, скорость
    res.gnss_t_meas    = g_t_meas(gidx);      % момент измерения
    res.gnss_time_age  = g_age(gidx);         % возраст измерения [с]
    res.gnss_mismatch  = g_mismatch(gidx,:);  % r_ant(t_meas) - r_ant(t_nav)
    res.gnss_offset    = gnss_dt_off;         % заданный сдвиг [с]
    res.init_pos_err   = init_pos_err;        % применённая начальная ошибка
    res.init_vel_err   = init_vel_err;
    res.gnss_n_shift   = n_shift;             % сдвиг в отсчётах сетки

    res.nav    = nav;
    res.P      = kf.P;

    % Окна аутэйджа, ОБРЕЗАННЫЕ по реальной длительности полёта
    % (иначе подсветка растягивает ось X далеко за конец траектории)
    res.outage = cfg.outage;
    res.outage(:,1) = max(cfg.outage(:,1), truth.t(1));
    res.outage(:,2) = min(cfg.outage(:,2), truth.t(end));

    % Ошибка в локальной ENU (для КВО)
    res.dr_enu   = res.dr * truth.Cen;

    % --- Ошибка в NED (внешний интерфейс; core остаётся ECEF) ---
    % Величины хранятся построчно, поэтому умножаем справа на NED -> ECEF.
    % Обратная совместимость: если truth сгенерирована старой версией и поля
    % C_e_ned0 нет, выводим матрицу из Cen перестановкой столбцов
    % (ned = [N E -U]); численно это тождественно прямому построению.
    if isfield(truth, 'C_e_ned0')
        C_e_ned0 = truth.C_e_ned0;
    else
        C_e_ned0 = [truth.Cen(:,2), truth.Cen(:,1), -truth.Cen(:,3)];
    end
    res.C_e_ned0 = C_e_ned0;
    res.dr_ned   = res.dr * C_e_ned0;
    res.dv_ned   = res.dv * C_e_ned0;
    res.err_horiz_final = hypot(res.dr_enu(end,1), res.dr_enu(end,2));
    res.err_3d_final    = norm(res.dr(end,:));
end


function ok = gnss_available(tk, outage)
%GNSS_AVAILABLE  Проверка доступности GNSS: false внутри окон аутэйджа.
    ok = true;
    for i = 1:size(outage,1)
        if tk >= outage(i,1) && tk < outage(i,2)
            ok = false;
            return;
        end
    end
end


function k = get_cal(prof, name)
%GET_CAL  Коэффициент остатка после предпусковой калибровки.
%   Раздельные cal_factor_gyro / cal_factor_accel имеют приоритет над общим
%   cal_factor; при отсутствии обоих принимается 1 (калибровки нет).
%   Логика повторяет imu_draw_errors, чтобы P0 строилась по тем же данным,
%   что и сама истина.
    if isfield(prof, name)
        k = prof.(name);
    elseif isfield(prof, 'cal_factor')
        k = prof.cal_factor;
    else
        k = 1.0;
    end
end
