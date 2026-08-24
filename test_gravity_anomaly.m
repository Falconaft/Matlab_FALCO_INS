%% TEST_GRAVITY_ANOMALY  Детерминированная проверка аномалии поля (B2)
%
%  ЦЕЛЬ: подтвердить, что известная постоянная аномалия dg вносится в
%  механизацию с правильным ЗНАКОМ и правильным МАСШТАБОМ.
%
%  ОЖИДАЕМОЕ ПОВЕДЕНИЕ. При dg = g_true - g_model ошибка ускорения
%  навигатора относительно истины равна -dg, поэтому на участке БЕЗ
%  коррекции, начиная с нулевых начальных ошибок:
%
%       dv(t) ~ -dg * t
%       dr(t) ~ -0.5 * dg * t^2
%
%  где dr = nav - truth (та же конвенция, что в res.dr).
%
%  ПОЧЕМУ ИНТЕРВАЛ КОРОТКИЙ. Простая аналитическая оценка верна лишь пока
%  не накопились эффекты второго порядка: Кориолис (-2[w_ie x]dv), градиент
%  гравитации по смещению позиции, и изменение ориентации. Их вклад растёт
%  быстрее, чем главный член, поэтому проверка ведётся на первых секундах
%  полёта, где отношение сигнал/искажение максимально.
%
%  УСЛОВИЯ ТЕСТА:
%    - все ошибки IMU нулевые (идеальный датчик);
%    - начальные ошибки нулевые;
%    - GNSS ОТКЛЮЧЁН на всём полёте (окно аутэйджа на всю траекторию);
%    - dg задан ЯВНО (не разыгрывается), чтобы сравнение было точным.
%
%  Требует falco_step7.mat.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг G.1: КОНФИГУРАЦИЯ
% =====================================================================
cfg = falco_config();
cfg.seed       = 4242;
cfg.f_diag     = cfg.fnav;          % диагностика на КАЖДОМ такте
cfg.diag_decim = 1;

t_end = truth.t(end);

% GNSS отключён полностью: контур разомкнут, ошибка накапливается свободно
cfg.outage = [ truth.t(1) - 1.0,  t_end + 1.0 ];

% Нулевые начальные ошибки
cfg.align_err   = zeros(3,1);
cfg.align_sigma = zeros(3,1);
cfg.P0_att      = [1e-3; 1e-3; 1e-3];
cfg.P0_pos      = [5.0; 5.0; 5.0];
cfg.P0_vel      = [0.1; 0.1; 0.1];

% =====================================================================
% Шаг G.2: ИДЕАЛЬНЫЙ ДАТЧИК
% =====================================================================
prof_ideal = zero_imu_errors(imu_profile_pulse40_updated());

rng_ideal = RandStream('mt19937ar','Seed', cfg.seed);
e_ideal   = imu_draw_errors(prof_ideal, rng_ideal);
imu_meas  = add_imu_errors(imu, prof_ideal, e_ideal, rng_ideal);

% =====================================================================
% Шаг G.3: ЗАДАННАЯ АНОМАЛИЯ (детерминированная, не розыгрыш)
% =====================================================================
g0       = 9.80665;
RAD2SEC  = 206265;

defl_E_sec = 10.0;      % уклонение отвеса на восток [угл.сек]
defl_N_sec = -6.0;      % уклонение отвеса на север  [угл.сек]
anom_U_mgal = 40.0;     % аномалия модуля g          [мГал]

dg_enu = [ g0 * defl_E_sec/RAD2SEC;      % Восток
           g0 * defl_N_sec/RAD2SEC;      % Север
           anom_U_mgal * 1e-5 ];         % Вверх
dg_ecef = truth.Cen * dg_enu;

cfg.dg_model = dg_ecef;

fprintf('=== TEST GRAVITY ANOMALY (B2) ===\n');
fprintf('Заданная аномалия (ENU):\n');
fprintf('  уклонение отвеса   : E %+.2f", N %+.2f"  (полное %.2f")\n', ...
        defl_E_sec, defl_N_sec, hypot(defl_E_sec, defl_N_sec));
fprintf('  аномалия вертикали : %+.1f мГал\n', anom_U_mgal);
fprintf('  dg в ENU  : [%+.4e %+.4e %+.4e] м/с²\n', dg_enu);
fprintf('  dg в ECEF : [%+.4e %+.4e %+.4e] м/с²\n', dg_ecef);
fprintf('  |dg| горизонтальная = %.4e м/с²\n', hypot(dg_enu(1), dg_enu(2)));
fprintf('  |dg| полная         = %.4e м/с²\n\n', norm(dg_ecef));

% =====================================================================
% Шаг G.4: ПРОГОН БЕЗ GNSS
% =====================================================================
t_run = tic;
res = eskf_run(imu_meas, truth, c, prof_ideal, cfg);
fprintf('Прогон завершён за %.1f с, обновлений GNSS: %d [%s]\n\n', ...
        toc(t_run), res.n_gnss, ...
        ternary_str(res.n_gnss == 0, 'OK — контур разомкнут', 'ОШИБКА'));

% =====================================================================
% Шаг G.5: СРАВНЕНИЕ С АНАЛИТИКОЙ НА КОРОТКОМ ИНТЕРВАЛЕ
% =====================================================================
% Аналитика:  dv = -dg*t,  dr = -0.5*dg*t^2   (dr = nav - truth)
t_check = [1.0, 2.0, 3.0, 5.0, 8.0, 12.0];

fprintf('%s\n', repmat('=',1,74));
fprintf('  СРАВНЕНИЕ С АНАЛИТИКОЙ  (dr = nav - truth, ожидаем -0.5*dg*t^2)\n');
fprintf('%s\n', repmat('=',1,74));
fprintf('  %5s | %11s %11s %8s | %11s %11s %8s\n', ...
        't, с', '|dv| факт', '|dv| теор', 'отн.', '|dr| факт', '|dr| теор', 'отн.');
fprintf('  %s\n', repmat('-',1,70));

for tq = t_check
    if tq > res.t(end), continue; end
    [~, k] = min(abs(res.t - tq));
    tk = res.t(k);

    dv_act = res.dv(k,:)';
    dr_act = res.dr(k,:)';
    dv_th  = -dg_ecef * tk;
    dr_th  = -0.5 * dg_ecef * tk^2;

    fprintf('  %5.1f | %11.3e %11.3e %8.4f | %11.3e %11.3e %8.4f\n', ...
            tk, norm(dv_act), norm(dv_th), norm(dv_act)/max(norm(dv_th),eps), ...
                norm(dr_act), norm(dr_th), norm(dr_act)/max(norm(dr_th),eps));
end
fprintf('%s\n', repmat('=',1,74));

% --- Покомпонентная сверка в ECEF на опорный момент ---
t_ref = 5.0;
[~, kr] = min(abs(res.t - t_ref));
tr = res.t(kr);
dv_th = -dg_ecef * tr;
dr_th = -0.5 * dg_ecef * tr^2;

fprintf('\nПОКОМПОНЕНТНО в ECEF при t = %.2f с:\n', tr);
fprintf('  dv факт : [%+.5e %+.5e %+.5e] м/с\n', res.dv(kr,:));
fprintf('  dv теор : [%+.5e %+.5e %+.5e] м/с\n', dv_th);
fprintf('  dr факт : [%+.5e %+.5e %+.5e] м\n',   res.dr(kr,:));
fprintf('  dr теор : [%+.5e %+.5e %+.5e] м\n',   dr_th);

err_dv = norm(res.dv(kr,:)' - dv_th) / max(norm(dv_th), eps);
err_dr = norm(res.dr(kr,:)' - dr_th) / max(norm(dr_th), eps);
fprintf('\n  относительное расхождение dv : %.4f  [%s]\n', err_dv, ...
        pass_str(err_dv < 0.05));
fprintf('  относительное расхождение dr : %.4f  [%s]\n', err_dr, ...
        pass_str(err_dr < 0.05));

% --- Проверка ЗНАКА покомпонентно ---
sign_ok = all(sign(res.dr(kr,:)') == sign(dr_th) | abs(dr_th) < 1e-12);
fprintf('  знак совпадает по всем осям  : %s\n', pass_str(sign_ok));

% =====================================================================
% Шаг G.6: КОНТРОЛЬ dg = 0 (регрессия)
% =====================================================================
fprintf('\n%s\n', repmat('=',1,74));
fprintf('  КОНТРОЛЬ: тот же прогон при dg = 0\n');
fprintf('%s\n', repmat('=',1,74));
cfg0 = cfg;
cfg0.dg_model = zeros(3,1);
res0 = eskf_run(imu_meas, truth, c, prof_ideal, cfg0);

[~, k0] = min(abs(res0.t - t_ref));
fprintf('  при t = %.2f с:  |dr| = %.4e м,  |dv| = %.4e м/с\n', ...
        t_ref, norm(res0.dr(k0,:)), norm(res0.dv(k0,:)));
fprintf('  terminal      :  |dr| = %.4e м,  |dv| = %.4e м/с\n', ...
        norm(res0.dr(end,:)), norm(res0.dv(end,:)));
fprintf('  (это численный пол механизации; должен совпасть с TEST A)\n');

fprintf('\n  Отношение с аномалией / без неё при t = %.2f с: %.1f раз\n', ...
        t_ref, norm(res.dr(kr,:))/max(norm(res0.dr(k0,:)), eps));

save('falco_test_gravity.mat', 'res', 'res0', 'dg_ecef', 'dg_enu', 'cfg');
fprintf('\nсохранено: falco_test_gravity.mat\n');


% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% =========================================================================

function prof_out = zero_imu_errors(prof_in)
%ZERO_IMU_ERRORS  Копия профиля со всеми обнулёнными источниками ошибок.
%   Модель ошибок НЕ переписывается: обнуляются только СКО в КОПИИ структуры.
    prof_out = prof_in;

    prof_out.gyro.turnon_bias_sigma  = 0;
    prof_out.gyro.inrun_bias_sigma   = 0;
    prof_out.gyro.arw                = 0;
    prof_out.gyro.scale_factor_sigma = 0;
    prof_out.gyro.misalign_axes      = 0;
    prof_out.gyro.misalign_ortho     = 0;
    prof_out.gyro.vrc                = 0;
    prof_out.gyro.g_sens_cal_sigma   = 0;

    prof_out.accel.turnon_bias_sigma  = 0;
    prof_out.accel.inrun_bias_sigma   = 0;
    prof_out.accel.vrw                = 0;
    prof_out.accel.scale_factor_sigma = 0;
    prof_out.accel.misalign_axes      = 0;
    prof_out.accel.misalign_ortho     = 0;
    prof_out.accel.vrc                = 0;

    prof_out.vib_grms = 0;
end


function s = pass_str(ok)
%PASS_STR  Метка результата проверки.
    if ok, s = 'OK'; else, s = 'FAIL'; end
end


function s = ternary_str(cond, a, b)
%TERNARY_STR  Выбор строки по условию.
    if cond, s = a; else, s = b; end
end
