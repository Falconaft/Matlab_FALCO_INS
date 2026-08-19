%% MAIN_STEP8  Шаг 8: один прогон ESKF на РЕАЛИСТИЧНОЙ траектории (Шаг 7)
%
%  Отличия от main_step5 (старая траектория):
%    - загружает falco_step7.mat (реалистичная траектория)
%    - окна аутэйджа вычисляются из p.coast_duration, а не заданы вручную
%    - удельная сила НЕ ноль на коасте -> ошибка ориентации ПРОЕЦИРУЕТСЯ
%      в позицию, поэтому главным источником становится точность ориентации,
%      а не bias акселерометра (оценка: 1 мрад -> 17.6 м за 60 с при |f|=1g)
%
%  Требует falco_step7.mat (запусти main_step7 один раз, генерация ~95 с).

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг 8.1: КОНФИГУРАЦИЯ
% =====================================================================
cfg.fnav       = 400;              % частота навигатора [Гц]
cfg.f_gnss     = 10;                % частота GNSS [Гц]
cfg.sigma_pos  = 5*ones(3,1);    % СКО шума позиции GNSS [м]
cfg.sigma_vel  = 0.1*ones(3,1);   % СКО шума скорости GNSS [м/с]
cfg.corrtime   = 1000;             % время корреляции Гаусса-Маркова [с]
cfg.lever      = [0.1; 0.0; -0.1]; % плечо антенны GNSS от IMU в body [м]
cfg.seed       = 3000;
cfg.diag_decim = 25;

% --- Окна аутэйджа: вычисляются из параметров траектории ---
t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;
cfg.outage = [ 0.0            p.outage_boost_end;   % аутэйдж на разгоне
               t_coast_start  t_end + 1.0 ];        % терминальный коаст

% --- Начальная выставка ---
cfg.align_err = [0.5e-3; 0.5e-3; 4e-3];   % ФАКТИЧЕСКАЯ ошибка [рад]
cfg.P0_att    = 1.3 * abs(cfg.align_err);   % что фильтр думает о ней
cfg.P0_pos    = [5.0; 5.0; 5.0];
cfg.P0_vel    = [0.1; 0.1; 0.1];

% =====================================================================
% Шаг 8.2: ПРОФИЛЬ IMU И ЗАШУМЛЁННЫЕ ИЗМЕРЕНИЯ
% =====================================================================
prof = imu_profile_pulse40_updated;

rng_imu  = RandStream('mt19937ar','Seed', cfg.seed);
e_imu    = imu_draw_errors(prof, rng_imu);
imu_meas = add_imu_errors(imu, prof, e_imu, rng_imu);

g0 = 9.80665;
fprintf('=== ШАГ 8: ESKF НА РЕАЛИСТИЧНОЙ ТРАЕКТОРИИ ===\n');
fprintf('Траектория: полёт %.1f с, апогей %.2f км, дальность %.2f км\n', ...
        t_end, max(truth.alt)/1000, hypot(truth.enu(end,1),truth.enu(end,2))/1000);
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

% =====================================================================
% Шаг 8.3: ПРОГОН ФИЛЬТРА
% =====================================================================
tic;
res = eskf_run(imu_meas, truth, c, prof, cfg);
fprintf('Прогон завершён за %.1f с\n\n', toc);

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
fprintf('\nГоризонтальная ошибка на ударе (КВО-метрика): %.2f м\n', res.err_horiz_final);
fprintf('Полная 3D-ошибка на ударе:                   %.2f м\n', res.err_3d_final);

% % --- Ключевая диагностика: ориентация на входе в коаст ---
% [~, k_coast] = min(abs(res.t - t_coast_start));
% dpsi_coast = norm(res.dpsi(k_coast,:))*1e3;
% fprintf('\n--- КРИТИЧНО ДЛЯ НОВОЙ ТРАЕКТОРИИ ---\n');
% fprintf('Ошибка ориентации на входе в коаст: %.3f мрад\n', dpsi_coast);
% fprintf('Оценка её вклада в позицию за коаст: %.2f м\n', ...
%         0.5*1.0*g0*dpsi_coast*1e-3*p.coast_duration^2);
% fprintf('(при |f|~1g; требование для 10 м: < 0.57 мрад)\n');

fprintf('\nОценка bias акс.: %s мг\n', mat2str(res.nav.ba'/g0*1e3, 4));
fprintf('Истинный bias акс.: %s мг\n', mat2str(e_imu.accel_turnon_bias'/g0*1e3, 4));
fprintf('Оценено: %.1f%%\n', ...
        100*(1 - norm(e_imu.accel_turnon_bias - res.nav.ba)/norm(e_imu.accel_turnon_bias)));

% =====================================================================
% Шаг 8.5: ГРАФИКИ
% =====================================================================
plot_eskf_results(res, truth, e_imu.accel_turnon_bias, e_imu.gyro_turnon_bias);

save('falco_step8.mat', 'res', 'cfg', 'prof');
fprintf('\nсохранено: falco_step8.mat\n');
