clear; clc;

load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

cfg  = falco_config();
prof = imu_profile_pulse40_updated();

N_mc      = 6;
base_seed = 7000;

% ============================================================
% КОНФИГУРАЦИЯ КАК В РЕАЛЬНОМ СЦЕНАРИИ
% ============================================================

cfg.sigma_pos = 5*ones(3,1);
cfg.sigma_vel = 0.1*ones(3,1);

cfg.f_diag     = 1;
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);

cfg.cep_target = 10;

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;

cfg.outage = [ ...
    0.0,           p.outage_boost_end;
    t_coast_start, t_end + 1.0 ];

cfg.align_sigma = [0.5e-3; 0.5e-3; 4.0e-3];
cfg.P0_att      = 1.3 * cfg.align_sigma;

cfg.init_pos_sigma = 1.0*ones(3,1);
cfg.init_vel_sigma = 0.03*ones(3,1);

cfg.P0_pos = cfg.init_pos_sigma;
cfg.P0_vel = cfg.init_vel_sigma;

% Аномалия гравитации
cfg.defl_vert_sigma = 10 / 206265;
cfg.grav_anom_sigma = 50e-5;

% Рассинхронизация GNSS/INS
cfg.gnss_time_offset_enable = true;
cfg.gnss_time_offset        = 2e-3;

% ============================================================
% ПРОФИЛЬ IMU
% ============================================================

prof.cal_factor_gyro  = 0.003;
prof.cal_factor_accel = 0.25;

% ============================================================
% SERIAL
% ============================================================

cfg.use_parallel = false;

mc_serial = run_montecarlo_eskf2( ...
    prof, prof, imu, truth, c, cfg, N_mc, base_seed);

% ============================================================
% PARALLEL
% ============================================================

cfg.use_parallel     = true;
cfg.parallel_workers = 6;

mc_parallel = run_montecarlo_eskf2( ...
    prof, prof, imu, truth, c, cfg, N_mc, base_seed);

% ============================================================
% СРАВНЕНИЕ
% ============================================================

d_horiz = max(abs( ...
    mc_serial.err_horiz - mc_parallel.err_horiz));

d_3d = max(abs( ...
    mc_serial.err_3d - mc_parallel.err_3d));

d_enu = max(abs( ...
    mc_serial.err_enu(:) - mc_parallel.err_enu(:)));

fprintf('\n=== SERIAL vs PARALLEL ===\n');
fprintf('max |delta err_horiz| = %.16g\n', d_horiz);
fprintf('max |delta err_3d|    = %.16g\n', d_3d);
fprintf('max |delta err_enu|   = %.16g\n', d_enu);

if d_horiz == 0 && d_3d == 0 && d_enu == 0
    fprintf('РЕЗУЛЬТАТ: ПОЛНОЕ СОВПАДЕНИЕ\n');
else
    fprintf('РЕЗУЛЬТАТ: есть различия\n');
end