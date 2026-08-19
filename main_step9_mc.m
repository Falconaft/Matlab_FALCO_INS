%% MAIN_STEP9_MC  Шаг 9: Монте-Карло на РЕАЛИСТИЧНОЙ траектории
%
%  Новый базовый КВО для Pulse-40 на траектории с ненулевой удельной силой.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг 9.1: ПАРАМЕТРЫ КАМПАНИИ
% =====================================================================
N_mc       = 25;         %  20, потом 200-250
base_seed  = 3001;
cep_target = 10;

% =====================================================================
% Шаг 9.2: КОНФИГУРАЦИЯ
% =====================================================================
cfg.fnav       = 250;
cfg.f_gnss     = 10;
cfg.sigma_pos  = 5*ones(3,1);
cfg.sigma_vel  = 0.1*ones(3,1);
cfg.corrtime   = 1000;
cfg.lever      = [0.1; 0.0; -0.1];
% Частота записи диагностики задаётся ЯВНО в герцах, прореживание
% вычисляется из неё. Раньше стояло число 250, привязанное к f_nav=250
% (250 тактов = 1 с); при смене частоты навигатора смысл менялся бы молча.
cfg.f_diag     = 1;                                   % [Гц]
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);
cfg.cep_target = cep_target;

% Аномалия гравитационного поля (навигатор её не знает).
% Уклонение отвеса — ГЛАВНОЕ для КВО (горизонтальная ошибка):
% 5 угл.сек равнина, 15 холмистая местность, 30 горы.
cfg.defl_vert_sigma = 10 / 206265;      % [рад] СКО уклонения отвеса
cfg.grav_anom_sigma = 50e-5;            % [м/с²] СКО аномалии (50 мГал)

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;
cfg.outage = [ 0.0            p.outage_boost_end;
               t_coast_start  t_end + 1.0 ];

% Выставка: СКО (разыгрывается на каждой реализации)
cfg.align_sigma = [0.5e-3; 0.5e-3; 4.0e-3];
cfg.P0_att      = 1.3 * cfg.align_sigma;
cfg.P0_pos      = [5.0; 5.0; 5.0];
cfg.P0_vel      = [0.1; 0.1; 0.1];

prof = imu_profile_pulse40_updated();

% =====================================================================
% ПРОВЕРКА ВРЕМЕННОЙ СЕТКИ (A1)
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

fprintf('=== ШАГ 9: МОНТЕ-КАРЛО НА РЕАЛИСТИЧНОЙ ТРАЕКТОРИИ ===\n');
fprintf('Кандидат: %s\n', prof.name);
fprintf('Траектория: %.1f с, коаст %.1f с, окно GNSS %.1f с\n', ...
        t_end, p.coast_duration, t_coast_start - p.outage_boost_end);
fprintf('Выставка (СКО): курс %.1f мрад\n', cfg.align_sigma(3)*1e3);
fprintf('Реализаций: %d\n\n', N_mc);

% =====================================================================
% Шаг 9.3: ПРОГОН
% =====================================================================
mc = run_montecarlo_eskf(prof, imu, truth, c, cfg, N_mc, base_seed);

% =====================================================================
% Шаг 9.4: РЕЗУЛЬТАТ
% =====================================================================
fprintf('\n=== РЕЗУЛЬТАТ: %s (реалистичная траектория) ===\n', mc.prof_name);
fprintf('%s\n', repmat('=',1,54));
fprintf('  КВО (медиана)               : %8.2f м\n', mc.cep);
fprintf('  Среднее                     : %8.2f м\n', mc.mean_horiz);
fprintf('  СКО                         : %8.2f м\n', mc.std_horiz);
fprintf('  R95                         : %8.2f м\n', mc.r95);
fprintf('  Худший прогон               : %8.2f м\n', mc.max_horiz);
fprintf('  Уложилось в %2.0f м            : %7.1f %%\n', cep_target, 100*mc.frac_pass);
fprintf('  Расходившихся               : %8d\n', mc.n_fail);
fprintf('%s\n', repmat('-',1,54));
fprintf('  Оценка bias акс. (медиана)  : %7.1f %%\n', 100*mc.ba_frac_est);
fprintf('  Систематика (ENU)           : [%+.2f %+.2f %+.2f] м\n', mc.bias_enu);
fprintf('%s\n', repmat('=',1,54));

fprintf('\nСРАВНЕНИЕ СО СТАРОЙ ТРАЕКТОРИЕЙ:\n');
fprintf('  старая (|f|=0 на коасте):  КВО ~1.4-1.9 м\n');
fprintf('  новая  (|f|=0.13..1.75g):  КВО  %.2f м\n', mc.cep);
if mc.cep > 3
    fprintf('  -> рост подтверждает: прежняя оценка была ЗАВЫШЕНА\n');
    fprintf('     из-за ненаблюдаемости ориентации на «стерильной» баллистике\n');
else
    fprintf('  -> КВО сопоставим: лучшая наблюдаемость компенсирует\n');
    fprintf('     проекцию ошибки ориентации в позицию\n');
end

plot_montecarlo_cep(mc, cep_target);
save('falco_step9_mc.mat', 'mc', 'cfg', 'prof', 'N_mc', 'base_seed');
fprintf('\nсохранено: falco_step9_mc.mat\n');
