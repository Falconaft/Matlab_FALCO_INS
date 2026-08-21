%% TEST_NAV_RATE_REGRESSION  Чистое сравнение f_nav = 250 и 400 Гц
%
%  ЦЕЛЬ:
%    Сравнить влияние ТОЛЬКО частоты навигатора на один и тот же
%    стохастический сценарий.
%
%  Ключевой принцип:
%    - truth одна и та же;
%    - профиль IMU один и тот же;
%    - одна и та же реализация turn-on bias / scale / misalignment / GM / ARW / VRW;
%    - одна и та же реализация GNSS-шума;
%    - одинаковые начальные ошибки, P0, outage, lever arm;
%    - меняется ТОЛЬКО cfg.fnav и связанный с ним diag_decim.
%
%  Поэтому это apples-to-apples regression test для 250 -> 400 Гц.
%
%  ТРЕБОВАНИЕ:
%    falco_step7.mat должен быть уже перегенерирован на truth = 2000 Гц
%    (truth.dt = 0.0005 с).
%
%  Ожидаем:
%    2000/250 = 8  — целое
%    2000/400 = 5  — целое
%    250/10   = 25 — целое
%    400/10   = 40 — целое

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг R.1: ОБЩАЯ КОНФИГУРАЦИЯ
% =====================================================================
cfg_common.f_gnss     = 10;
cfg_common.sigma_pos  = 5*ones(3,1);
cfg_common.sigma_vel  = 0.1*ones(3,1);
cfg_common.corrtime   = 1000;
cfg_common.lever      = [0.1; 0.0; -0.1];
cfg_common.seed       = 3000;
cfg_common.f_diag     = 10;

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;

cfg_common.outage = [ 0.0            p.outage_boost_end;
                      t_coast_start  t_end + 1.0 ];

cfg_common.align_err = [0.5e-3; 0.5e-3; 4e-3];
cfg_common.P0_att    = 1.3 * abs(cfg_common.align_err);
cfg_common.P0_pos    = [5.0; 5.0; 5.0];
cfg_common.P0_vel    = [0.1; 0.1; 0.1];

prof = imu_profile_pulse40_updated();

% =====================================================================
% Шаг R.2: ПРОВЕРКА СЕТКИ
% =====================================================================
f_truth = 1/truth.dt;

check_rate(f_truth, 250, cfg_common.f_gnss);
check_rate(f_truth, 400, cfg_common.f_gnss);

fprintf('\n=== TEST_NAV_RATE_REGRESSION ===\n');
fprintf('truth = %.0f Гц, GNSS = %.0f Гц\n', f_truth, cfg_common.f_gnss);
fprintf('Сравнение выполняется на ОДНОЙ и той же реализации IMU/GNSS.\n\n');

% =====================================================================
% Шаг R.3: ОДНА СТОХАСТИЧЕСКАЯ РЕАЛИЗАЦИЯ IMU
% =====================================================================
rng_imu  = RandStream('mt19937ar','Seed', cfg_common.seed);
e_imu    = imu_draw_errors(prof, rng_imu);
imu_meas = add_imu_errors(imu, prof, e_imu, rng_imu);

g0 = 9.80665;
fprintf('--- ОБЩАЯ РЕАЛИЗАЦИЯ IMU ---\n');
fprintf('seed                : %d\n', cfg_common.seed);
fprintf('turn-on bias gyro   : %.4f °/ч\n', ...
        norm(e_imu.gyro_turnon_bias)/(pi/180)*3600);
fprintf('turn-on bias accel  : %.4f мг\n', ...
        norm(e_imu.accel_turnon_bias)/g0*1e3);

% =====================================================================
% Шаг R.4: ПРОГОН 250 Гц
% =====================================================================
cfg250 = cfg_common;
cfg250.fnav       = 250;
cfg250.diag_decim = round(cfg250.fnav / cfg250.f_diag);

fprintf('\n=========================================================\n');
fprintf('  RUN A: f_nav = 250 Гц\n');
fprintf('=========================================================\n');
t250 = tic;
res250 = eskf_run(imu_meas, truth, c, prof, cfg250);
runtime250 = toc(t250);

m250 = collect_metrics(res250, e_imu, runtime250, g0);

% =====================================================================
% Шаг R.5: ПРОГОН 400 Гц
% =====================================================================
cfg400 = cfg_common;
cfg400.fnav       = 400;
cfg400.diag_decim = round(cfg400.fnav / cfg400.f_diag);

fprintf('\n=========================================================\n');
fprintf('  RUN B: f_nav = 400 Гц\n');
fprintf('=========================================================\n');
t400 = tic;
res400 = eskf_run(imu_meas, truth, c, prof, cfg400);
runtime400 = toc(t400);

m400 = collect_metrics(res400, e_imu, runtime400, g0);

% =====================================================================
% Шаг R.6: СВОДНОЕ СРАВНЕНИЕ
% =====================================================================
fprintf('\n==========================================================================\n');
fprintf('  СРАВНЕНИЕ 250 vs 400 Гц  (ОДНА И ТА ЖЕ реализация IMU/GNSS)\n');
fprintf('==========================================================================\n');
fprintf('  %-28s %14s %14s %14s\n', 'Метрика', '250 Гц', '400 Гц', '400-250');
fprintf('  %s\n', repmat('-',1,76));

print_row('terminal horiz, м', m250.err_horiz, m400.err_horiz);
print_row('terminal 3D, м',    m250.err_3d,    m400.err_3d);
print_row('terminal |dr|, м',  m250.term_dr,   m400.term_dr);
print_row('terminal |dv|, м/с',m250.term_dv,   m400.term_dv);
print_row('terminal |dpsi|, мрад', m250.term_dpsi_mrad, m400.term_dpsi_mrad);
print_row('max |dr|, м',       m250.max_dr,    m400.max_dr);
print_row('max |dv|, м/с',     m250.max_dv,    m400.max_dv);
print_row('max |dpsi|, мрад',  m250.max_dpsi_mrad, m400.max_dpsi_mrad);

fprintf('  %-28s %14.2f %14.2f %14.2f\n', ...
        'bias accel estimated, %', ...
        m250.ba_frac_est*100, m400.ba_frac_est*100, ...
        (m400.ba_frac_est - m250.ba_frac_est)*100);

fprintf('  %-28s %14d %14d %14d\n', ...
        'GNSS updates', m250.n_gnss, m400.n_gnss, ...
        m400.n_gnss - m250.n_gnss);

fprintf('  %-28s %14.2f %14.2f %14.2f\n', ...
        'runtime, с', m250.runtime, m400.runtime, ...
        m400.runtime - m250.runtime);

fprintf('==========================================================================\n');

fprintf('\n--- ОТНОШЕНИЯ 400/250 ---\n');
fprintf('  horiz error ratio : %.4f\n', safe_ratio(m400.err_horiz, m250.err_horiz));
fprintf('  3D error ratio    : %.4f\n', safe_ratio(m400.err_3d, m250.err_3d));
fprintf('  |dr| ratio        : %.4f\n', safe_ratio(m400.term_dr, m250.term_dr));
fprintf('  |dv| ratio        : %.4f\n', safe_ratio(m400.term_dv, m250.term_dv));
fprintf('  |dpsi| ratio      : %.4f\n', safe_ratio(m400.term_dpsi_mrad, m250.term_dpsi_mrad));
fprintf('  runtime ratio     : %.4f\n', safe_ratio(m400.runtime, m250.runtime));

% =====================================================================
% Шаг R.7: СОХРАНЕНИЕ
% =====================================================================
regression.f_truth = f_truth;
regression.prof     = prof;
regression.e_imu    = e_imu;
regression.cfg250   = cfg250;
regression.cfg400   = cfg400;
regression.res250   = res250;
regression.res400   = res400;
regression.m250     = m250;
regression.m400     = m400;

save('falco_nav_rate_regression.mat', 'regression');
fprintf('\nсохранено: falco_nav_rate_regression.mat\n');


% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% =========================================================================

function check_rate(f_truth, f_nav, f_gnss)
%CHECK_RATE  Fail-fast проверка кратности сеток.

    ratio_truth_nav = f_truth / f_nav;
    ratio_nav_gnss  = f_nav / f_gnss;

    assert(abs(ratio_truth_nav - round(ratio_truth_nav)) < 1e-12, ...
        'test_nav_rate_regression:badTruthNavRatio', ...
        'f_truth/f_nav = %.12f не целое для f_nav = %.3f Гц.', ...
        ratio_truth_nav, f_nav);

    assert(abs(ratio_nav_gnss - round(ratio_nav_gnss)) < 1e-12, ...
        'test_nav_rate_regression:badNavGnssRatio', ...
        'f_nav/f_gnss = %.12f не целое для f_nav = %.3f Гц.', ...
        ratio_nav_gnss, f_nav);

    fprintf('Проверка сетки: truth/nav = %.0f/%.0f = %.0f, nav/GNSS = %.0f/%.0f = %.0f [OK]\n', ...
        f_truth, f_nav, ratio_truth_nav, f_nav, f_gnss, ratio_nav_gnss);
end


function m = collect_metrics(res, e_imu, runtime, g0)
%COLLECT_METRICS  Метрики одного прогона.

    dr_n   = vecnorm(res.dr,   2, 2);
    dv_n   = vecnorm(res.dv,   2, 2);
    dpsi_n = vecnorm(res.dpsi, 2, 2) * 1e3;  % мрад

    m.err_horiz      = res.err_horiz_final;
    m.err_3d         = res.err_3d_final;
    m.term_dr        = dr_n(end);
    m.term_dv        = dv_n(end);
    m.term_dpsi_mrad = dpsi_n(end);

    m.max_dr         = max(dr_n);
    m.max_dv         = max(dv_n);
    m.max_dpsi_mrad  = max(dpsi_n);

    m.n_gnss         = res.n_gnss;
    m.runtime        = runtime;

    if norm(e_imu.accel_turnon_bias) > 0
        m.ba_frac_est = 1 - ...
            norm(e_imu.accel_turnon_bias - res.nav.ba) / ...
            norm(e_imu.accel_turnon_bias);
    else
        m.ba_frac_est = NaN;
    end

    m.ba_est_mg  = res.nav.ba / g0 * 1e3;
    m.ba_true_mg = e_imu.accel_turnon_bias / g0 * 1e3;
end


function print_row(name, a, b)
%PRINT_ROW  Одна строка таблицы сравнения.
    fprintf('  %-28s %14.6e %14.6e %14.6e\n', name, a, b, b-a);
end


function r = safe_ratio(a, b)
%SAFE_RATIO  Отношение с защитой от деления на ноль.
    if abs(b) < eps
        r = NaN;
    else
        r = a / b;
    end
end
