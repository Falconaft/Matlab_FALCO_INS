function p = imu_profile_pulse40()
%IMU_PROFILE_PULSE40  Профиль ошибок IMU SBG Systems Pulse-40 (OEM).
%   Значения из OEM User Manual v1 (Performance specifications §5).
%   Спецификации даны как typical (1σ) по всему диапазону -40..+85°C.
%   Все значения переведены в СИ.

    deg = pi/180;
    g0  = 9.80665;
    p.name = 'SBG Pulse-40';

    % =====================================================================
    % ГИРОСКОП (§5.2)
    % =====================================================================
    % In-run bias instability (Allan, const T) = 0.8 °/ч
    p.gyro.inrun_bias_sigma     = 0.8 * deg/3600;      % [рад/с]
    p.gyro.inrun_bias_corr_time = 1000;                % [с]  <- ОЦЕНКА
    % Long-term bias repeatability = 250 °/ч (turn-on)
    p.gyro.turnon_bias_sigma    = 250 * deg/3600;      % [рад/с]
    % Angular Random Walk (Allan, const T) = 0.08 °/√ч
    p.gyro.arw                  = 0.08 * deg/sqrt(3600); % [рад/√с]
    % Scale Factor error = 1000 ppm (для ±490°/с range)
    p.gyro.scale_factor_sigma   = 1000e-6;             % [-]
    % Acceleration sensitivity = 10 °/ч/g (datasheet, tested over ±1g)
    p.gyro.g_sensitivity        = 10 * deg/3600 / g0;  % [(рад/с)/(м/с²)]
    % Misalignment = 0.34 mrad
    p.gyro.misalignment         = 0.34e-3;             % [рад]

    % =====================================================================
    % АКСЕЛЕРОМЕТР ±40g (§5.1)
    % =====================================================================
    % In-run bias instability (Allan, const T) = 6 µg
    p.accel.inrun_bias_sigma     = 6e-6 * g0;          % [м/с²]
    p.accel.inrun_bias_corr_time = 1000;               % [с]  <- ОЦЕНКА
    % Long-term bias repeatability = 1250 µg (turn-on)
    % (сноска: экспортная стандартная версия имеет bias > 1250 µg)
    p.accel.turnon_bias_sigma    = 1250e-6 * g0;       % [м/с²]
    % Velocity Random Walk = 0.02 м/с/√ч
    p.accel.vrw                  = 0.02 / sqrt(3600);  % [м/с/√с]
    % Scale Factor error = 300 ppm
    p.accel.scale_factor_sigma   = 300e-6;             % [-]
    % Misalignment = 0.34 mrad
    p.accel.misalignment         = 0.34e-3;            % [рад]

    % =====================================================================
    % ЭКСПЛУАТАЦИЯ
    % =====================================================================
    p.cal_factor      = 1.0;
    p.accel.range_g   = 40;                            % OK: штатно ±40g
    p.max_sample_rate = 6600;                          % gyro sampling 6.6 kHz
end