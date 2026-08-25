function imu_noisy = add_imu_errors(imu_ideal, profile, e, rng_stream)
%ADD_IMU_ERRORS  Наложение полной модели ошибок IMU на идеальные показания.
%
%   Модель по Groves ("Principles of GNSS...", eq.4.16-4.17) с добавленной
%   динамикой in-run bias (Гаусс-Марков 1-го порядка):
%
%     ω_изм = (I + diag(s_g))·ω_ист + b_turnon + b_inrun(t) + G_g_true·f_ист + шум
%     f_изм = (I + diag(s_a))·f_ист + b_turnon + b_inrun(t) + шум
%
%   G_g_true = G_g_nom .* (1 + d) — реальная g-чувствительность экземпляра,
%   где d — ошибка калибровки из imu_draw_errors. Механизация компенсирует
%   её ПАСПОРТНЫМ G_g_nom, поэтому в навигацию проходит только остаток.
%
%   Вход:
%     imu_ideal  - идеальные показания (.t, .fb, .wib_b, .Ceb)
%     profile    - профиль IMU
%     e          - реализация постоянных ошибок (imu_draw_errors)
%     rng_stream - поток случайных чисел
%   Выход:
%     imu_noisy  - зашумлённые .fb, .wib_b (та же сетка)

    t   = imu_ideal.t;
    fb  = imu_ideal.fb;
    wib = imu_ideal.wib_b;
    N   = numel(t);
    dt  = t(2) - t(1);

    % =====================================================================
    % Шаг 1: ПАРАМЕТРЫ ГАУССА-МАРКОВА для in-run bias
    % =====================================================================
    % b[k] = beta*b[k-1] + sigma*sqrt(1-beta^2)*eta,  eta~N(0,1)
    % Такой выбор возбуждающего шума даёт стационарное СКО = sigma.
    Tg   = profile.gyro.inrun_bias_corr_time;
    Ta   = profile.accel.inrun_bias_corr_time;
    beta_g = exp(-dt/Tg);
    beta_a = exp(-dt/Ta);
    qd_g = profile.gyro.inrun_bias_sigma  * sqrt(1 - beta_g^2);
    qd_a = profile.accel.inrun_bias_sigma * sqrt(1 - beta_a^2);

    % Начальное значение из стационарного распределения
    b_inrun_g = profile.gyro.inrun_bias_sigma  * randn(rng_stream, 3, 1);
    b_inrun_a = profile.accel.inrun_bias_sigma * randn(rng_stream, 3, 1);

    % =====================================================================
    % Шаг 2: МАТРИЦЫ SCALE FACTOR и g-SENSITIVITY
    % =====================================================================
    % ПОЛНЫЕ матрицы ошибок: масштаб + поворот триады + неортогональность.
    % Внедиагональные члены (перекос осей) фильтр НЕ оценивает — в состоянии
    % только диагональный масштаб. Это делает их чистой неучтённой ошибкой,
    % вносящей паразитную поперечную составляющую: сила вдоль X даёт
    % показание по Y и Z.
    Ma = eye(3) + imu_error_matrix(e.accel_scale, e.accel_misalign, e.accel_ortho);
    Mg = eye(3) + imu_error_matrix(e.gyro_scale,  e.gyro_misalign,  e.gyro_ortho);

    % Паспортная (номинальная) g-чувствительность — её же использует
    % механизация для компенсации.
    Gg_nom = profile.gyro.g_sensitivity * eye(3);

    % Реальная g-чувствительность экземпляра: номинал с ошибкой калибровки.
    Gg_true = Gg_nom * diag(1 + e.g_sens_cal_err);

    % Коэффициенты белого шума: root-PSD -> СКО дискретного шума за такт
    arw_std = profile.gyro.arw  / sqrt(dt);
    vrw_std = profile.accel.vrw / sqrt(dt);

    % =====================================================================
    % ВИБРАЦИОННЫЙ РЕКТИФИКАЦИОННЫЙ КОЭФФИЦИЕНТ (VRC)
    % =====================================================================
    % Под вибрацией нелинейность датчика выпрямляет знакопеременный сигнал
    % и создаёт ЛОЖНОЕ ПОСТОЯННОЕ смещение, пропорциональное КВАДРАТУ уровня
    % вибрации: b_vrc = VRC * g_rms^2.
    %
    % Тонкость: вибрация есть только на работающем двигателе. Фильтр за это
    % время успевает "выучить" ложное смещение, а после отсечки вибрация
    % пропадает — остаётся НЕВЕРНАЯ оценка bias, которая портит коаст.
    % Это реальный эффект, а не артефакт модели.
    vrc_a = get_field(profile.accel, 'vrc', 0);   % [(м/с²)/g²]
    vrc_g = get_field(profile.gyro,  'vrc', 0);   % [(рад/с)/g²]
    vib_grms   = get_field(profile, 'vib_grms',   0);    % уровень вибрации [g rms]
    vib_t_end  = get_field(profile, 'vib_t_end',  0);    % до какого времени [с]

    % Смещение = СКО (модуль из datasheet) * случайный множитель экземпляра.
    % Множитель разыгран в imu_draw_errors и постоянен в течение полёта;
    % между реализациями Монте-Карло он меняется, поэтому систематики
    % по ансамблю не возникает (в отличие от прежней версии с ones(3,1)).
    vrc_dir_a = get_vec(e, 'accel_vrc_dir');
    vrc_dir_g = get_vec(e, 'gyro_vrc_dir');
    b_vrc_a = vrc_a * vib_grms^2 * vrc_dir_a;
    b_vrc_g = vrc_g * vib_grms^2 * vrc_dir_g;

    fb_n  = zeros(N,3);
    wib_n = zeros(N,3);

    % ИСТОРИЯ ПОЛНОГО ИСТИННОГО СМЕЩЕНИЯ (только для диагностики).
    % Модель НЕ меняется: сохраняются величины, которые и так вычисляются
    % в цикле. Нужны потому, что фильтр оценивает ПОЛНОЕ смещение
    % (turn-on + in-run), и сравнивать его качество с одним лишь turn-on
    % некорректно: у гироскопа Pulse-40 in-run bias (0.5 °/ч) РАВЕН
    % калиброванному turn-on (150*0.0033 = 0.495 °/ч), поэтому такое
    % сравнение давало отрицательные сотни процентов.
    bg_total = zeros(N,3);      % полное истинное смещение гироскопа [рад/с]
    ba_total = zeros(N,3);      % полное истинное смещение акселерометра [м/с²]

    % =====================================================================
    % Шаг 3: ОСНОВНОЙ ЦИКЛ — наложение ошибок такт за тактом
    % =====================================================================
    for k = 1:N
        f_true = fb(k,:)';
        w_true = wib(k,:)';

        % --- (а) Эволюция in-run bias (Гаусс-Марков) ---
        b_inrun_g = beta_g*b_inrun_g + qd_g*randn(rng_stream,3,1);
        b_inrun_a = beta_a*b_inrun_a + qd_a*randn(rng_stream,3,1);

        % --- (б) Белый шум (ARW/VRW) ---
        noise_g = arw_std * randn(rng_stream,3,1);
        noise_a = vrw_std * randn(rng_stream,3,1);

        % --- Вибрация действует только при работающем двигателе ---
        if t(k) <= vib_t_end
            vrc_bias_a = b_vrc_a;
            vrc_bias_g = b_vrc_g;
        else
            vrc_bias_a = zeros(3,1);
            vrc_bias_g = zeros(3,1);
        end

        % --- (в) Акселерометр ---
        f_meas = Ma*f_true + e.accel_turnon_bias + b_inrun_a + vrc_bias_a + noise_a;

        % --- (г) Гироскоп (+ РЕАЛЬНАЯ g-чувствительность экземпляра) ---
        w_meas = Mg*w_true + e.gyro_turnon_bias + b_inrun_g ...
                 + Gg_true*f_true + vrc_bias_g + noise_g;

        fb_n(k,:)  = f_meas';
        wib_n(k,:) = w_meas';

        % Полное истинное смещение на этом такте (диагностика).
        % VRC включён: фильтр не отличает его от bias и оценивает вместе.
        bg_total(k,:) = (e.gyro_turnon_bias  + b_inrun_g + vrc_bias_g)';
        ba_total(k,:) = (e.accel_turnon_bias + b_inrun_a + vrc_bias_a)';
    end

    % =====================================================================
    % Шаг 4: УПАКОВКА
    % =====================================================================
    imu_noisy.t      = t;
    imu_noisy.fb     = fb_n;
    imu_noisy.wib_b  = wib_n;
    imu_noisy.Ceb    = imu_ideal.Ceb;
    imu_noisy.Gg_nom = Gg_nom;      % передаём в механизацию для компенсации

    % СЛУЖЕБНАЯ ИСТИННАЯ КИНЕМАТИКА (не показание датчика!).
    % w_eb_b — истинная угловая скорость тела относительно ECEF, сформированная
    % в attitude_program. Пробрасывается БЕЗ ИЗМЕНЕНИЙ и БЕЗ ЗАШУМЛЕНИЯ,
    % поскольку нужна для синтеза ИСТИННОЙ скорости GNSS-антенны:
    %     v_ant_true = v_IMU_true + C_true*(w_eb_true x lever)
    % Загрязнять её ошибками гироскопа нельзя: это часть ИСТИНЫ, а не
    % измерение. Механизация и фильтр это поле не используют.
    imu_noisy.w_eb_b = imu_ideal.w_eb_b;

    % История полного истинного смещения — только для диагностики качества
    % оценки bias. Механизация и фильтр эти поля НЕ используют.
    imu_noisy.bg_total = bg_total;
    imu_noisy.ba_total = ba_total;
end


function v = get_field(s, name, default)
%GET_FIELD  Чтение поля структуры со значением по умолчанию.
%   Обеспечивает обратную совместимость со старыми профилями.
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end


function v = get_vec(s, name)
%GET_VEC  Чтение 3x1 поля структуры; нули, если поля нет.
%   Обратная совместимость со старыми реализациями ошибок.
    if isfield(s, name)
        v = s.(name);
    else
        v = zeros(3,1);
    end
end