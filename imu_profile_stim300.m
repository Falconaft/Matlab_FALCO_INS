function p = imu_profile_stim300()
%IMU_PROFILE_STIM300  Профиль ошибок IMU Sensonor/Safran STIM300 (версия 30g).
%   Значения из datasheet TS1524 rev.26. Акселерометр - конфигурация ±30g
%   (Table 6-6), выбрана потому, что пиковая ОСЕВАЯ нагрузка на разгоне ~16g
%   требует диапазона с запасом (10g-версия НАСЫЩАЕТСЯ -> катастрофа, проверено).
%   Гироскопные параметры общие для всех версий акселерометра.
%
%   Все значения переведены в СИ. Пометка "<- ОЦЕНКА" там, где datasheet не даёт
%   параметр напрямую.
%
%   ЕДИНИЦЫ datasheet -> СИ:
%     °/h -> рад/с: *(pi/180)/3600 ; °/√h -> рад/√с: *(pi/180)/sqrt(3600)
%     mg -> м/с²: *1e-3*g0 ; m/s/√h -> м/с/√с: /sqrt(3600)
%     °/h/g -> (рад/с)/(м/с²): *(pi/180)/3600/g0

    deg = pi/180;
    g0  = 9.80665;
    p.name = 'STIM300 (30g)';

    % =====================================================================
    % ГИРОСКОП (общий для всех версий акселерометра)
    % =====================================================================
    % Bias instability (Allan) ~0.3 °/h (по Allan-графику datasheet)
    p.gyro.inrun_bias_sigma     = 0.3 * deg/3600;      % [рад/с]  <- по Allan-графику
    % Время корреляции in-run bias: datasheet не даёт -> оценка
    p.gyro.inrun_bias_corr_time = 1000;                % [с]  <- ОЦЕНКА
    % Bias error over temperature (run-run) = 4 °/h -> turn-on repeatability
    p.gyro.turnon_bias_sigma    = 4 * deg/3600;        % [рад/с]
    % Angular Random Walk = 0.15 °/√h (datasheet, Allan @25°C)
    p.gyro.arw                  = 0.15 * deg/sqrt(3600); % [рад/√с]
    % Scale factor accuracy = ±5 % (datasheet; трактуем как 1σ, консервативно)
    p.gyro.scale_factor_sigma   = 0.05;                % [-]
    % g-sensitivity (bias, с g-компенсацией) = 1 °/h/g (datasheet Linear Accel Effect)
    p.gyro.g_sensitivity        = 1 * deg/3600 / g0;   % [(рад/с)/(м/с²)]
    % Misalignment = 1 mrad (datasheet)
    p.gyro.misalignment         = 1e-3;                % [рад]

    % =====================================================================
    % АКСЕЛЕРОМЕТР — версия ±30g (Table 6-6 datasheet)
    % =====================================================================
    % Bias Instability (Allan @25°C) = 0.12 mg
    p.accel.inrun_bias_sigma     = 0.12e-3 * g0;       % [м/с²]
    % Время корреляции in-run bias: datasheet не даёт -> оценка
    p.accel.inrun_bias_corr_time = 1000;               % [с]  <- ОЦЕНКА
    % Turn-on bias: берём "Bias error over temperature ±6 mg rms" как 1σ
    % (консервативнее, чем switch on/off ±2.3 mg границы). См. примечание ниже.
    p.accel.turnon_bias_sigma    = 6e-3 * g0;          % [м/с²]
    % Velocity Random Walk (Allan @25°C) = 0.21 m/s/√h
    p.accel.vrw                  = 0.21 / sqrt(3600);  % [м/с/√с]
    % Scale Factor 1 year stability = 300 ppm (nom)
    p.accel.scale_factor_sigma   = 300e-6;             % [-]
    % Misalignment = 1 mrad (Orthogonality 0.6 mrad - отдельно, здесь не используется)
    p.accel.misalignment         = 1e-3;               % [рад]

    % ПРИМЕЧАНИЕ по turn-on акселерометра:
    %   datasheet даёт ДВА параметра:
    %     - "Bias switch on/off repeatability" ±2.3 mg (границы, ~3σ -> 1σ≈0.77 mg)
    %     - "Bias error over temperature" ±6 mg rms (это уже 1σ)
    %   Взят второй (6 mg) как более полный и консервативный. Для оптимистичного
    %   сценария (термостабилизация) можно снизить до 0.77 mg.

    % =====================================================================
    % ПАРАМЕТРЫ КАЛИБРОВКИ / ЭКСПЛУАТАЦИИ
    % =====================================================================
    % Доля turn-on bias, остающаяся после стартовой калибровки (0..1).
    p.cal_factor      = 1.0;
    % Диапазон акселерометра [g]: ±30g покрывает пик 16g с запасом ×1.9.
    p.accel.range_g   = 30;
    % Максимальная частота выборки [Гц].
    p.max_sample_rate = 2000;
end