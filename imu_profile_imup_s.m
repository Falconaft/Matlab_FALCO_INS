function p = imu_profile_imup_s()
%IMU_PROFILE_IMUP_S  Профиль ошибок IMU Inertial Labs IMU-P Tactical S.
%   Значения из datasheet Rev 1.4 (Apr 2024). АКСЕЛЕРОМЕТР — конфигурация ±40g
%   (не ±8g!), выбрана из-за пиковой осевой нагрузки ~16g на разгоне.
%   Все значения переведены в СИ. "<- ОЦЕНКА" где datasheet не даёт параметр.

    deg = pi/180;
    g0  = 9.80665;
    p.name = 'IMU-P Tactical S';

    % =====================================================================
    % ГИРОСКОП
    % =====================================================================
    % Bias in-run stability (Allan) = 0.5 °/ч
    p.gyro.inrun_bias_sigma     = 0.5 * deg/3600;      % [рад/с]
    p.gyro.inrun_bias_corr_time = 1000;                % [с]  <- ОЦЕНКА
    % Bias repeatability (turn-on to turn-on) = 15 °/ч  (НЕ over-temp 35!)
    p.gyro.turnon_bias_sigma    = 15 * deg/3600;       % [рад/с]
    % Angular Random Walk = 0.06 °/√ч
    p.gyro.arw                  = 0.06 * deg/sqrt(3600); % [рад/√с]
    % SF accuracy (over temperature) = 300 ppm
    p.gyro.scale_factor_sigma   = 300e-6;              % [-]
    % g-sensitivity: НЕ указана -> оценка 1 °/ч/g (уровень tactical, как STIM300)
    p.gyro.g_sensitivity        = 1 * deg/3600 / g0;   % [(рад/с)/(м/с²)]  <- ОЦЕНКА
    % Axis misalignment = 0.15 mrad
    p.gyro.misalignment         = 0.15e-3;             % [рад]

    % =====================================================================
    % АКСЕЛЕРОМЕТР — версия ±40g (правая колонка datasheet)
    % =====================================================================
    % Bias in-run stability (Allan) = 0.03 mg
    p.accel.inrun_bias_sigma     = 0.03e-3 * g0;       % [м/с²]
    p.accel.inrun_bias_corr_time = 1000;               % [с]  <- ОЦЕНКА
    % Turn-on: bias one-year repeatability = 1.5 mg (см. примечание про over-temp)
    p.accel.turnon_bias_sigma    = 1.5e-3 * g0;        % [м/с²]
    % Velocity Random Walk = 0.045 м/с/√ч
    p.accel.vrw                  = 0.045 / sqrt(3600); % [м/с/√с]
    % SF accuracy (over temperature) = 500 ppm
    p.accel.scale_factor_sigma   = 500e-6;             % [-]
    % Axis misalignment = 0.15 mrad
    p.accel.misalignment         = 0.15e-3;            % [рад]

    % ПРИМЕЧАНИЕ turn-on акселерометра ±40g:
    %   datasheet: "bias one-year repeatability" 1.5 mg, "bias instability over
    %   temperature RMS" 1.2 mg. Взят 1.5 mg (repeatability). Для консервативного
    %   теплового сценария можно взять 1.2 mg over-temp.

    % =====================================================================
    % ЭКСПЛУАТАЦИЯ
    % =====================================================================
    p.cal_factor      = 1.0;
    p.accel.range_g   = 40;                            % OK: покрывает пик 16g (×2.5)
    p.max_sample_rate = 4000;                          % data update rate 4000 Hz
end