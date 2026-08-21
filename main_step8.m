%% MAIN_STEP8  Шаг 8: один прогон ESKF на РЕАЛИСТИЧНОЙ траектории (Шаг 7)
%
%  Отличия от main_step5 (старая траектория):
%    - загружает falco_step7.mat (реалистичная траектория)
%    - окна аутэйджа вычисляются из p.coast_duration, а не заданы вручную
%    - удельная сила НЕ ноль на коасте -> ошибка ориентации ПРОЕЦИРУЕТСЯ
%      в позицию
%
%  РЕЖИМ 400 Гц (Шаг B). Сетка истины повышена до 2000 Гц в traj_params.m,
%  поэтому 2000/400 = 5 — целое, и фактическая частота навигатора точно
%  равна заявленной. До этого при сетке 1000 Гц отношение 1000/400 = 2.5
%  округлялось до 3, и реальная частота составляла 333.3 Гц.
%
%  ТРЕБУЕТСЯ ПЕРЕГЕНЕРАЦИЯ falco_step7.mat запуском main_step7 после
%  изменения p.dt, иначе truth.dt останется прежним и make_increments
%  прервёт выполнение по assert.
%
%  Опорные значения BASELINE_250 (f_truth=1000, f_nav=250) для сравнения:
%     горизонтальная ошибка 2.193 м, 3D 2.311 м, runtime 11.36 с,
%     эпох навигатора 52143.
%
%  Требует falco_step7.mat (запусти main_step7 один раз).

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг 8.1: КОНФИГУРАЦИЯ
% =====================================================================
% Базовая конфигурация; ниже переопределяется только специфичное для сценария.
cfg = falco_config();

% --- Специфика одиночного прогона Шага 8 ---
cfg.seed = 3000;
% sigma_pos, sigma_vel, P0_pos совпадают с базовыми (5 м, 0.1 м/с, 5 м),
% поэтому не переопределяются.

% Частота записи диагностики задаётся ЯВНО в герцах, а прореживание
% вычисляется из неё. Раньше здесь стояло число 25, привязанное к 250 Гц
% (25 тактов = 0.1 с); при смене частоты навигатора такая запись молча
% меняла бы смысл. Теперь семантика сохраняется при любой fnav.
cfg.f_diag     = 10;                                  % [Гц]
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);   % пересчёт после смены f_diag

% --- Окна аутэйджа: вычисляются из параметров траектории ---
t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;
cfg.outage = [ 0.0            p.outage_boost_end;   % аутэйдж на разгоне
               t_coast_start  t_end + 1.0 ];        % терминальный коаст

% --- Начальная выставка ---
cfg.align_err = [1.5e-3; 1.5e-3; 15e-3];    % ФАКТИЧЕСКАЯ ошибка [рад]
cfg.P0_att    = 1.3 * abs(cfg.align_err);  % что фильтр думает о ней
cfg.P0_pos    = [5.0; 5.0; 5.0];
cfg.P0_vel    = [0.1; 0.1; 0.1];

% =====================================================================
% Шаг 8.1b: ПРОВЕРКА ВРЕМЕННОЙ СЕТКИ (A1)
% =====================================================================
f_truth   = 1/truth.dt;
ratio_raw = f_truth / cfg.fnav;
step_chk  = round(ratio_raw);
dt_nav    = step_chk * truth.dt;
fnav_real = 1/dt_nav;

fprintf('=== ПРОВЕРКА ВРЕМЕННОЙ СЕТКИ (A1) ===\n');
fprintf('  запрошенная f_nav   : %10.4f Гц\n', cfg.fnav);
fprintf('  f_truth             : %10.4f Гц\n', f_truth);
fprintf('  ratio f_truth/f_nav : %10.6f\n',    ratio_raw);
fprintf('  step (отсчётов/такт): %10d\n',      step_chk);
fprintf('  dt_nav              : %10.6f с\n',  dt_nav);
fprintf('  ФАКТИЧЕСКАЯ f_nav   : %10.4f Гц\n', fnav_real);
fprintf('  ratio f_nav/f_gnss  : %10.4f\n',    cfg.fnav/cfg.f_gnss);
if abs(ratio_raw - step_chk) > 1e-12
    fprintf('  СТАТУС: НЕКРАТНЫЕ ЧАСТОТЫ — реально %.3f Гц вместо %.3f Гц!\n', ...
            fnav_real, cfg.fnav);
else
    fprintf('  СТАТУС: OK — частоты кратны, f_nav точная\n');
end
fprintf('\n');

% =====================================================================
% Шаг 8.2: ПРОФИЛЬ IMU И ЗАШУМЛЁННЫЕ ИЗМЕРЕНИЯ
% =====================================================================
prof = imu_profile_pulse40_updated();

rng_imu  = RandStream('mt19937ar','Seed', cfg.seed);
e_imu    = imu_draw_errors(prof, rng_imu);
imu_meas = add_imu_errors(imu, prof, e_imu, rng_imu);

g0 = 9.80665;
fprintf('=== ШАГ 8: ESKF НА РЕАЛИСТИЧНОЙ ТРАЕКТОРИИ (f_nav = %d Гц) ===\n', cfg.fnav);
fprintf('Траектория: полёт %.1f с, апогей %.2f км, дальность %.2f км\n', ...
        t_end, max(truth.alt)/1000, hypot(truth.enu(end,1),truth.enu(end,2))/1000);
fprintf('Отсчётов истины: %d, эпох навигатора: ~%d\n', ...
        numel(truth.t), floor((numel(truth.t)-1)/step_chk));
fprintf('Профиль IMU: %s\n', prof.name);
fprintf('Истинный turn-on bias: гиро %.1f °/ч, акс %.2f мг\n', ...
        norm(e_imu.gyro_turnon_bias)/(pi/180)*3600, ...
        norm(e_imu.accel_turnon_bias)/g0*1e3);
fprintf('Выставка: крен/тангаж %.1f мрад, курс %.1f мрад\n', ...
        cfg.align_err(1)*1e3, cfg.align_err(3)*1e3);
fprintf('GNSS: аутэйдж 0..%.1f с, активен %.1f..%.1f с (окно %.1f с)\n', ...
        p.outage_boost_end, p.outage_boost_end, t_coast_start, ...
        t_coast_start - p.outage_boost_end);
fprintf('Терминальный коаст: %.1f с\n\n', p.coast_duration);

% --- Характеристики траектории (для протокола BASELINE) ---
w_mag = vecnorm(imu.wib_b, 2, 2);
f_mag = vecnorm(imu.fb,    2, 2);
fprintf('--- ХАРАКТЕРИСТИКИ ТРАЕКТОРИИ ---\n');
fprintf('  апогей            : %.3f км\n',  max(truth.alt)/1000);
fprintf('  дальность         : %.3f км\n',  hypot(truth.enu(end,1),truth.enu(end,2))/1000);
fprintf('  макс. скорость    : %.2f м/с\n', max(truth.Vmag));
fprintf('  макс. |f|         : %.4f g\n',   max(f_mag)/g0);
fprintf('  макс. |omega|     : %.4f °/с\n', rad2deg(max(w_mag)));
fprintf('  terminal truth r  : [%.3f %.3f %.3f] м\n', truth.R(end,:));
fprintf('  terminal truth v  : [%.4f %.4f %.4f] м/с\n', truth.V(end,:));
fprintf('\n');

% =====================================================================
% Шаг 8.3: ПРОГОН ФИЛЬТРА
% =====================================================================
t_run = tic;
res = eskf_run(imu_meas, truth, c, prof, cfg);
runtime = toc(t_run);
fprintf('Прогон завершён за %.2f с\n\n', runtime);

% =====================================================================
% Шаг 8.4: МЕТРИКИ
% =====================================================================
t_marks = [p.outage_boost_end, p.t_burn, t_coast_start-0.1, res.t(end)];
lbl     = {'конец аутэйджа разгона','конец работы двигателя', ...
           'перед потерей GNSS','удар'};
fprintf('%-26s %10s %10s %12s\n','Момент','|dr|, м','|dv|, м/с','|dpsi|, мрад');
fprintf('%s\n', repmat('-',1,62));
for i = 1:numel(t_marks)
    [~, k] = min(abs(res.t - t_marks(i)));
    fprintf('%-26s %10.2f %10.3f %12.3f\n', lbl{i}, ...
            norm(res.dr(k,:)), norm(res.dv(k,:)), norm(res.dpsi(k,:))*1e3);
end
fprintf('%s\n', repmat('-',1,62));
% --- Ошибка в NED (внешний интерфейс) ---
fprintf('\nОшибка позиции в NED на ударе: [%+.3f %+.3f %+.3f] м (N, E, D)\n', ...
        res.dr_ned(end,:));
fprintf('Горизонтальная ошибка на ударе (КВО-метрика): %.3f м\n', res.err_horiz_final);
fprintf('Полная 3D-ошибка на ударе:                   %.3f м\n', res.err_3d_final);

fprintf('\nОценка bias акс.: %s мг\n', mat2str(res.nav.ba'/g0*1e3, 4));
fprintf('Истинный bias акс.: %s мг\n', mat2str(e_imu.accel_turnon_bias'/g0*1e3, 4));
fprintf('Оценено: %.1f%%\n', ...
        100*(1 - norm(e_imu.accel_turnon_bias - res.nav.ba)/norm(e_imu.accel_turnon_bias)));

% =====================================================================
% Шаг 8.4b: ПРОТОКОЛ BASELINE_250 (сводка для сравнения с 400 Гц)
% =====================================================================
baseline.tag            = sprintf('RUN_%dHz', round(cfg.fnav));
baseline.f_truth        = f_truth;
baseline.f_nav_req      = cfg.fnav;
baseline.f_nav_real     = fnav_real;
baseline.dt_nav         = dt_nav;
baseline.step           = step_chk;
baseline.f_gnss         = cfg.f_gnss;
baseline.n_truth        = numel(truth.t);
baseline.n_nav_epochs   = floor((numel(truth.t)-1)/step_chk);
baseline.flight_time    = t_end;
baseline.apogee_km      = max(truth.alt)/1000;
baseline.range_km       = hypot(truth.enu(end,1),truth.enu(end,2))/1000;
baseline.v_max          = max(truth.Vmag);
baseline.f_max_g        = max(f_mag)/g0;
baseline.w_max_dps      = rad2deg(max(w_mag));
baseline.r_term_truth   = truth.R(end,:);
baseline.v_term_truth   = truth.V(end,:);
baseline.term_dr        = norm(res.dr(end,:));
baseline.term_dv        = norm(res.dv(end,:));
baseline.term_dpsi_mrad = norm(res.dpsi(end,:))*1e3;
baseline.err_horiz      = res.err_horiz_final;
baseline.err_3d         = res.err_3d_final;
baseline.runtime        = runtime;
baseline.seed           = cfg.seed;

fprintf('\n=========================================================\n');
fprintf('  ПРОТОКОЛ %s\n', baseline.tag);
fprintf('=========================================================\n');
fprintf('  f_truth / f_nav / f_gnss : %.0f / %.0f / %.0f Гц\n', ...
        baseline.f_truth, baseline.f_nav_real, baseline.f_gnss);
fprintf('  dt_nav / step            : %.6f с / %d\n', baseline.dt_nav, baseline.step);
fprintf('  отсчётов истины          : %d\n',   baseline.n_truth);
fprintf('  эпох навигатора          : %d\n',   baseline.n_nav_epochs);
fprintf('  время полёта             : %.3f с\n', baseline.flight_time);
fprintf('  апогей / дальность       : %.3f км / %.3f км\n', ...
        baseline.apogee_km, baseline.range_km);
fprintf('  макс V / |f| / |omega|   : %.2f м/с / %.4f g / %.4f °/с\n', ...
        baseline.v_max, baseline.f_max_g, baseline.w_max_dps);
fprintf('  terminal |dr|            : %.4f м\n',   baseline.term_dr);
fprintf('  terminal |dv|            : %.5f м/с\n', baseline.term_dv);
fprintf('  terminal |dpsi|          : %.4f мрад\n',baseline.term_dpsi_mrad);
fprintf('  горизонтальная ошибка    : %.4f м\n',   baseline.err_horiz);
fprintf('  3D-ошибка                : %.4f м\n',   baseline.err_3d);
fprintf('  runtime                  : %.2f с\n',   baseline.runtime);
fprintf('=========================================================\n');

% =====================================================================
% Шаг 8.5: ГРАФИКИ
% =====================================================================
plot_eskf_results(res, truth, e_imu.accel_turnon_bias, e_imu.gyro_turnon_bias);

% Сравнение истинной и навигационной траекторий (внешний NED-интерфейс)
plot_nav_comparison(truth, res, p, cfg, c);

save('falco_step8.mat', 'res', 'cfg', 'prof', 'baseline');
fprintf('\nсохранено: falco_step8.mat (включая структуру baseline)\n');
