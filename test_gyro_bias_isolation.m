%% TEST_GYRO_BIAS_ISOLATION
% Проверка: почему ESKF плохо оценивает gyro bias перед terminal coast.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

cfg = falco_config();

cfg.f_diag     = 1;
cfg.diag_decim = round(cfg.fnav/cfg.f_diag);
cfg.cep_target = 10;

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;

cfg.outage = [0.0            p.outage_boost_end;
              t_coast_start  t_end + 1.0];

cfg.align_sigma = [0.5e-3; 0.5e-3; 4.0e-3];
cfg.P0_att      = 1.3*cfg.align_sigma;
cfg.P0_pos      = [5;5;5];
cfg.P0_vel      = [0.1;0.1;0.1];

% B2/B4 оставляем как в baseline Step 9
cfg.defl_vert_sigma = 10/206265;
cfg.grav_anom_sigma = 50e-5;

cfg.gnss_time_offset_enable = true;
cfg.gnss_time_offset        = 2e-3;

prof_base = imu_profile_pulse40_updated();

N_mc      = 5;
base_seed = 3001;

%% --------------------------------------------------------------
% Сценарии
% prof_truth меняем, prof_filter ВСЕГДА остаётся prof_base.
%% --------------------------------------------------------------

names = {
    'BASELINE'
    'без gyro misalignment'
    'без gyro g-sens residual'
    'без gyro VRC'
    'только gyro bias + noise'
    'ЧИСТАЯ проверка gyro bias'
};

profiles = cell(6,1);

% 1. Всё включено
profiles{1} = prof_base;

% 2. Без misalignment/неортогональности
profiles{2} = prof_base;
profiles{2}.gyro.misalign_axes  = 0*prof_base.gyro.misalign_axes;
profiles{2}.gyro.misalign_ortho = 0*prof_base.gyro.misalign_ortho;

% 3. Без ошибки компенсации g-sensitivity
profiles{3} = prof_base;
profiles{3}.gyro.g_sens_cal_sigma = 0;

% 4. Без VRC гироскопа
profiles{4} = prof_base;
profiles{4}.gyro.vrc = 0;

% 5. Одновременно убираем неучтённые систематики gyro.
% Оставляем turn-on bias, in-run GM bias и ARW.
profiles{5} = prof_base;
profiles{5}.gyro.misalign_axes     = 0*prof_base.gyro.misalign_axes;
profiles{5}.gyro.misalign_ortho    = 0*prof_base.gyro.misalign_ortho;
profiles{5}.gyro.g_sens_cal_sigma  = 0;
profiles{5}.gyro.vrc               = 0;
profiles{5}.gyro.scale_factor_sigma = 0;

% 6. ЧИСТАЯ ПРОВЕРКА НАБЛЮДАЕМОСТИ GYRO BIAS
% Оставляем только:
%   - gyro turn-on bias
%   - gyro in-run GM bias
%   - gyro ARW
%   - GNSS measurement noise
%
% Убираем остальные ошибки, которые могут маскироваться под bg.

profiles{6} = prof_base;

% --- gyro: убрать всё кроме bias + ARW ---
profiles{6}.gyro.scale_factor_sigma = 0;

profiles{6}.gyro.misalign_axes  = 0 * prof_base.gyro.misalign_axes;
profiles{6}.gyro.misalign_ortho = 0 * prof_base.gyro.misalign_ortho;

profiles{6}.gyro.g_sens_cal_sigma = 0;
profiles{6}.gyro.vrc              = 0;

% --- accel: максимально очистить ---
profiles{6}.accel.turnon_bias_sigma = 0;
profiles{6}.accel.inrun_bias_sigma  = 0;
profiles{6}.accel.scale_factor_sigma = 0;

profiles{6}.accel.misalign_axes  = 0 * prof_base.accel.misalign_axes;
profiles{6}.accel.misalign_ortho = 0 * prof_base.accel.misalign_ortho;

profiles{6}.accel.vrc = 0;
profiles{6}.accel.vrw = 0;

%% --------------------------------------------------------------
% Прогоны
%% --------------------------------------------------------------

fprintf('\n=== GYRO BIAS ISOLATION TEST ===\n');
fprintf('N = %d, одинаковые seeds во всех сценариях\n\n', N_mc);

result = struct([]);

deg = pi/180;

for j = 1:numel(names)

    fprintf('\n[%d/%d] %s\n', j, numel(names), names{j});

    cfg_run = cfg;

% Для чистого теста gyro bias убираем внешние источники
if j == 6
    % Без B2
    cfg_run.defl_vert_sigma = 0;
    cfg_run.grav_anom_sigma = 0;

    % Без B4
    cfg_run.gnss_time_offset_enable = false;
    cfg_run.gnss_time_offset        = 0;

    % Без ошибки начальной выставки
    cfg_run.align_sigma = [0; 0; 0];

    % GNSS noise ОСТАВЛЯЕМ как в baseline
end

mc = run_montecarlo_eskf2( ...
    profiles{j}, ...     % truth errors
    prof_base, ...       % фильтр НЕ меняем
    imu, truth, c, cfg_run, N_mc, base_seed);

    true_bg = mc.bg_true_pre_med/deg*3600;
    resid   = mc.bg_resid_pre_med/deg*3600;

    result(j).name     = names{j};
    result(j).true_bg  = true_bg;
    result(j).resid_bg = resid;
    result(j).pre_psi  = mc.pre_dpsi_med;
    result(j).pre_dv   = mc.pre_dv_med;
    result(j).cep      = mc.cep;

    fprintf('  gyro true pre : %.3f deg/h\n', true_bg);
    fprintf('  gyro residual : %.3f deg/h\n', resid);
    fprintf('  pre |dpsi|    : %.3f mrad\n', mc.pre_dpsi_med);
    fprintf('  pre |dv|      : %.4f m/s\n', mc.pre_dv_med);
    fprintf('  CEP50         : %.3f m\n', mc.cep);
end

%% --------------------------------------------------------------
% Сводка
%% --------------------------------------------------------------

fprintf('\n\n===============================================================\n');
fprintf('%-28s %10s %10s %10s %10s\n', ...
    'Сценарий','true bg','resid bg','dpsi','CEP50');
fprintf('%-28s %10s %10s %10s %10s\n', ...
    '','deg/h','deg/h','mrad','m');
fprintf('---------------------------------------------------------------\n');

for j = 1:numel(result)
    fprintf('%-28s %10.3f %10.3f %10.3f %10.3f\n', ...
        result(j).name, ...
        result(j).true_bg, ...
        result(j).resid_bg, ...
        result(j).pre_psi, ...
        result(j).cep);
end

fprintf('===============================================================\n');