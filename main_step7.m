%% MAIN_STEP7  Шаг 7: генерация РЕАЛИСТИЧНОЙ траектории и синтез IMU
%
%  Профиль: разгон -> отсечка двигателя -> перелом программы тангажа ->
%           апогей -> баллистика -> квазибаллистический спуск -> терминал.
%
%  Ключевое отличие от Шага 1-2: удельная сила НЕ обнуляется после разгона.
%  На всех участках |f| = 0.13..1.8 g вместо машинного нуля, поэтому ошибка
%  ориентации и scale factor становятся НАБЛЮДАЕМЫМИ через GNSS.
%
%  Результат -> falco_step7.mat (совместим по интерфейсу с falco_step12.mat)

clear; clc; close all;

c = constants();
p = traj_params();

fprintf('=== ШАГ 7: РЕАЛИСТИЧНАЯ ТРАЕКТОРИЯ ===\n');
fprintf('Пусковая: %.4f°N %.4f°E, азимут %.0f°, возвышение %.0f°\n', ...
        p.lat0, p.lon0, p.azimuth, p.launch_el);
fprintf('Двигатель: тяга %.0f Н, время %.1f с, топливо %.0f кг\n', ...
        p.thrust, p.t_burn, p.m_fuel);
fprintf('Крен: амплитуда %.1f°, частота %.2f Гц\n\n', p.roll_amp_deg, p.roll_freq_hz);

% =====================================================================
% Шаг 7.A: ГЕНЕРАЦИЯ ТРАЕКТОРИИ
% =====================================================================
tic;
truth = gen_truth_v2(p, c);
fprintf('Траектория сгенерирована за %.1f с\n', toc);

N = numel(truth.t);
[apogee, i_ap] = max(truth.alt);
range_km = hypot(truth.enu(end,1), truth.enu(end,2))/1000;

fprintf('  Точек: %d, полёт %.2f с\n', N, truth.t(end));
fprintf('  Апогей: %.2f км на t=%.1f с\n', apogee/1000, truth.t(i_ap));
fprintf('  Дальность: %.2f км\n', range_km);
fprintf('  Макс. скорость: %.0f м/с (M=%.2f)\n', max(truth.Vmag), max(truth.Mach));
fprintf('  Конечная скорость: %.0f м/с, угол %.1f°\n\n', truth.Vmag(end), truth.gamma(end));

% =====================================================================
% Шаг 7.B: ПРОГРАММА ОРИЕНТАЦИИ И СИНТЕЗ IMU
% =====================================================================
tic;
imu = attitude_program(truth, p, c);
fprintf('IMU синтезирован за %.1f с\n', toc);

g0 = 9.80665;
f_mag = vecnorm(imu.fb, 2, 2)/g0;
w_mag = vecnorm(imu.wib_b, 2, 2);

fprintf('  Пиковая осевая перегрузка: %.2f g\n', max(abs(imu.fb(:,1)))/g0);
fprintf('  Пиковый модуль |f|: %.2f g\n', max(f_mag));
fprintf('  Макс. угловая скорость: %.2f °/с\n', rad2deg(max(w_mag)));
fprintf('  Суммарный поворот корпуса: %.1f°\n\n', rad2deg(sum(w_mag)*truth.dt));

% =====================================================================
% Шаг 7.C: ПРОВЕРКИ КОРРЕКТНОСТИ
% =====================================================================
fprintf('--- ПРОВЕРКИ ---\n');

% (1) Боковая удельная сила должна быть машинным нулём
% (движение в вертикальной плоскости, боковых сил нет)
fy_max = max(abs(imu.fb(:,2)))/g0;
fprintf('  |f_y| макс (должно быть ~0): %.2e g  [%s]\n', fy_max, ...
        ternary(fy_max < 1e-10, 'OK', 'ПРОВЕРИТЬ'));

% (2) Ортогональность матриц ориентации
orth_err = 0;
for k = 1:round(N/50):N
    orth_err = max(orth_err, norm(imu.Ceb(:,:,k)'*imu.Ceb(:,:,k) - eye(3)));
end
fprintf('  Ортогональность C_b_e: %.2e  [%s]\n', orth_err, ...
        ternary(orth_err < 1e-12, 'OK', 'ПРОВЕРИТЬ'));

% (3) Удельная сила НЕ обнуляется (главное отличие от старой модели)
i_after = truth.t > p.t_burn + 5;
fprintf('  |f| после разгона: %.3f .. %.3f g  [%s]\n', ...
        min(f_mag(i_after)), max(f_mag(i_after)), ...
        ternary(min(f_mag(i_after)) > 0.01, 'НАБЛЮДАЕМО', 'слишком мало'));

% (4) Гладкость: максимальная скорость изменения |f|
df = abs(diff(f_mag))./diff(truth.t);
fprintf('  Макс. d|f|/dt: %.1f g/с (при разрыве было бы >1000)\n\n', max(df));

% =====================================================================
% Шаг 7.D: УДЕЛЬНАЯ СИЛА ПО УЧАСТКАМ (наблюдаемость)
% =====================================================================
fprintf('--- УДЕЛЬНАЯ СИЛА ПО ФАЗАМ (определяет наблюдаемость) ---\n');
for ph = 1:4
    m = truth.phase == ph;
    if any(m)
        fprintf('  фаза %d: t=%6.1f..%6.1f с, |f| = %6.3f..%6.2f g, alpha=%+.1f°\n', ...
                ph, truth.t(find(m,1)), truth.t(find(m,1,'last')), ...
                min(f_mag(m)), max(f_mag(m)), truth.alpha(find(m,1,'last')));
    end
end

% =====================================================================
% Шаг 7.E: СЦЕНАРИЙ GNSS
% =====================================================================
t_coast_start = truth.t(end) - p.coast_duration;
fprintf('\n--- СЦЕНАРИЙ GNSS ---\n');
fprintf('  Аутэйдж на разгоне: 0.0 .. %.1f с\n', p.outage_boost_end);
fprintf('  GNSS активен: %.1f .. %.1f с (окно %.1f с)\n', ...
        p.outage_boost_end, t_coast_start, t_coast_start - p.outage_boost_end);
fprintf('  Терминальный коаст: %.1f .. %.1f с (%.1f с)\n', ...
        t_coast_start, truth.t(end), p.coast_duration);

% =====================================================================
% % Шаг 7.F: ГРАФИКИ
% % =====================================================================
% plot_trajectory_v2(truth, imu, p);

save('falco_step7.mat', 'truth', 'imu', 'c', 'p');
fprintf('\nсохранено: falco_step7.mat\n');


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
