function p = imu_profile_pulse80()
%IMU_PROFILE_PULSE80  Профиль ошибок IMU SBG Systems Pulse-80 (Box).
%   Значения из Box User Manual v2 (Performance specifications §5, PRELIMINARY).
%   Топовый датчик серии: гироскоп существенно лучше Pulse-40.
%   Спецификации typical (1σ), диапазон -40..+71°C. Все значения в СИ.

    deg = pi/180;
    g0  = 9.80665;
    p.name = 'SBG Pulse-80';

    % =====================================================================
    % ГИРОСКОП (§5.2)
    % =====================================================================
    % In-run bias instability (Allan, const T) = 0.1 °/ч
    p.gyro.inrun_bias_sigma     = 0.1 * deg/3600;      % [рад/с]
    p.gyro.inrun_bias_corr_time = 1000;                % [с]  <- ОЦЕНКА
    % Long-term bias repeatability = 20 °/ч (turn-on)
    p.gyro.turnon_bias_sigma    = 20 * deg/3600;       % [рад/с]
    % Angular Random Walk (Allan, const T) = 0.012 °/√ч
    p.gyro.arw                  = 0.012 * deg/sqrt(3600); % [рад/√с]
    % Scale Factor error = 150 ppm
    p.gyro.scale_factor_sigma   = 150e-6;              % [-]
    % g-sensitivity: в datasheet только VRC (0.08 °/h/g² - квадратичный, НЕ линейный).
    % Линейная g-sens не указана -> оценка как Pulse-40 (10 °/ч/g).
    p.gyro.g_sensitivity        = 10 * deg/3600 / g0;  % [(рад/с)/(м/с²)]  <- ОЦЕНКА
    % Misalignment = 0.34 mrad
    p.gyro.misalignment         = 0.34e-3;             % [рад]

    % =====================================================================
    % АКСЕЛЕРОМЕТР (§5.1)  — ПРИМЕЧАНИЕ: ±40g ТОЛЬКО ПО ЗАПРОСУ (export controlled)
    % =====================================================================
    % In-run bias instability (Allan, const T) = 6 µg
    p.accel.inrun_bias_sigma     = 6e-6 * g0;          % [м/с²]
    p.accel.inrun_bias_corr_time = 1000;               % [с]  <- ОЦЕНКА
    % Long-term bias repeatability = 1250 µg (turn-on)
    p.accel.turnon_bias_sigma    = 1250e-6 * g0;       % [м/с²]
    % Velocity Random Walk = 0.02 м/с/√ч
    p.accel.vrw                  = 0.02 / sqrt(3600);  % [м/с/√с]
    % Scale Factor error = 300 ppm
    p.accel.scale_factor_sigma   = 300e-6;             % [-]
    % Cross-axis sensitivity 0.01° -> misalignment ~0.17 mrad; берём 0.34 как Pulse
    p.accel.misalignment         = 0.34e-3;            % [рад]  <- ОЦЕНКА

    % =====================================================================
    % ЭКСПЛУАТАЦИЯ
    % =====================================================================
    p.cal_factor      = 1.0;
    % ВНИМАНИЕ: базовая версия НЕ ±40g. Диапазон ±40g доступен только по запросу
    % и export-controlled. Если базовый диапазон < 20g -> НАСЫЩЕНИЕ на разгоне.
    p.accel.range_g   = 40;                            % при условии заказа 40g версии
    p.accel.range_note = 'export-controlled, только по запросу';
    p.max_sample_rate = 4000;
end