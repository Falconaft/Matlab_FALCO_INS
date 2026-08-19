function p = imu_profile_pulse40_updated()
%IMU_PROFILE_PULSE40  Профиль ошибок IMU SBG Systems Pulse-40 Box.
%   Значения из Preliminary Specifications rev. 0.4 (§5).
%   Характеристики typical. Всё в СИ.

    deg = pi/180;
    g0  = 9.80665;
    p.name = 'SBG Pulse-40';

    % =====================================================================
    % ГИРОСКОП (§5.1)
    % =====================================================================
    p.gyro.inrun_bias_sigma     = 0.5 * deg/3600;       % in-run bias [рад/с]
    p.gyro.inrun_bias_corr_time = 1000;                 % [с] <- ОЦЕНКА
    p.gyro.turnon_bias_sigma    = 150 * deg/3600;       % bias repeatability
    p.gyro.arw                  = 0.08 * deg/sqrt(3600);% Angular Random Walk
    p.gyro.scale_factor_sigma   = 500e-6;               % Scale Factor error

    % g-чувствительность (Acceleration sensitivity, tested over ±1g)
    p.gyro.g_sensitivity        = 10 * deg/3600 / g0;   % [(рад/с)/(м/с²)]
    p.gyro.g_sens_cal_sigma     = 0.20;                 % 20%, допущение модели

    % Перекос осей
    p.gyro.misalign_axes        = 1.0e-3;               % Misalignment [рад]
    p.gyro.misalign_ortho       = 0.07e-3;              % Orthogonality [рад]

    % Вибрационный ректификационный коэффициент (VRC)
    % datasheet: 0.03 °/ч/g² RMS, 10g RMS, 20 Гц..2 кГц
    p.gyro.vrc                  = 0.03 * deg/3600;      % [(рад/с)/g²]

    % =====================================================================
    % АКСЕЛЕРОМЕТР ±40g (§5.2)
    % =====================================================================
    p.accel.inrun_bias_sigma     = 6e-6 * g0;           % in-run bias [м/с²]
    p.accel.inrun_bias_corr_time = 1000;                % [с] <- ОЦЕНКА
    p.accel.turnon_bias_sigma    = 1250e-6 * g0;        % bias repeatability
    p.accel.vrw                  = 0.02 / sqrt(3600);   % Velocity Random Walk
    p.accel.scale_factor_sigma   = 500e-6;              % Scale Factor error

    % Перекос осей
    p.accel.misalign_axes        = 1.0e-3;              % Misalignment [рад]
    p.accel.misalign_ortho       = 0.07e-3;             % Orthogonality [рад]

    % VRC: datasheet 0.03 mg/g² RMS
    p.accel.vrc                  = 0.03e-3 * g0;        % [(м/с²)/g²]

    % =====================================================================
    % ЭКСПЛУАТАЦИЯ И КАЛИБРОВКА
    % =====================================================================
    % Коэффициенты ниже — допущения модели, не параметры datasheet.
    p.cal_factor        = 1.0;       % общий fallback
    p.cal_factor_gyro   = 0.0033;    % ~0.5/150
    p.cal_factor_accel  = 0.25;      % предпусковая калибровка

    % Вибрация при работе двигателя — параметры сценария
    p.vib_grms   = 3.0;              % [g RMS]
    p.vib_t_end  = 8.0;              % [с]

    % Диапазоны и частоты
    p.accel.range_g             = 40;    % [g]
    p.gyro.range_dps            = 499;  % [°/с]
    p.gyro.internal_sample_rate = 6600;  % [Гц]
    p.accel.internal_sample_rate= 4000;  % [Гц]
    p.max_sample_rate           = 6600;  % [Гц]
end