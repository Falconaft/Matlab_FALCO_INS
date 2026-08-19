%% MAIN_STEP12  Falco NAV - Шаги 1 и 2: траектория и синтез ИНС
%
%  Шаг 1: генерация идеальной истинной траектории в ECEF (gen_truth).
%  Шаг 2: обратная механизация - синтез идеальных показаний ИНС (inverse_mech).
%
%  Два контрольных теста:
%    (A) Кажущаяся сила на баллистике ~ 0 (свободное падение -> акселерометр молчит).
%    (B) Угловая скорость ограничена и гладкая (нет вырождения программы ориентации).

clear; clc;

% =====================================================================
% КОНФИГУРАЦИЯ СЦЕНАРИЯ (Вариант 1: апогей 20 км, дальность 80 км)
% =====================================================================
c = constants();

cfg.lat0 = deg2rad(50.0);    % широта точки старта [рад]
cfg.lon0 = deg2rad(30.0);    % долгота точки старта [рад]
cfg.h0   = 0.0;              % высота старта [м]
cfg.t_b  = 10.0;             % длительность разгона [с]
cfg.v_up = 576.7;            % вертикальная (Up) скорость в конце разгона [м/с]
cfg.v_h  = 625.6;            % горизонтальная (East) скорость в конце разгона [м/с]
cfg.dt   = 1e-3;             % шаг интегрирования траектории [с] (1 кГц)
cfg.T    = 140.0;            % максимальное время (траектория обрежется на ударе) [с]

% =====================================================================
% ШАГ 1: ТРАЕКТОРИЯ
% =====================================================================
truth = gen_truth(cfg, c);

alt = (truth.R - truth.r0') * truth.Cen(:,3);   % высота над стартом (проекция на Up)
fprintf('=== ШАГ 1: ИСТИННАЯ ТРАЕКТОРИЯ ===\n');
fprintf(' Апогей           = %.4f км\n', max(alt)/1000);
fprintf(' Время до удара    = %.2f с   (N=%d точек)\n', truth.t_impact, numel(truth.t));
fprintf(' |V| конца разгона = %.3f м/с\n', norm(truth.V(truth.nb+1,:)));

% =====================================================================
% ШАГ 2: ОБРАТНАЯ МЕХАНИЗАЦИЯ (синтез ИНС)
% =====================================================================
imu = inverse_mech(truth, c);

f_norm = vecnorm(imu.fb, 2, 2);        % модуль кажущейся силы [м/с^2]
w_norm = vecnorm(imu.wib_b, 2, 2);     % модуль угловой скорости [рад/с]
nb = truth.nb;

fprintf('\n=== ШАГ 2: ОБРАТНАЯ МЕХАНИЗАЦИЯ ===\n');
fprintf(' ТЕСТ A |f^b| баллистика: max=%.3e м/с^2  (должно быть ~0)\n', max(f_norm(nb+6:end)));
fprintf('        |f^b| разгон:     max=%.2f  mean=%.2f м/с^2\n', max(f_norm(1:nb)), mean(f_norm(1:nb)));
fprintf(' ТЕСТ B |w_ib^b| разгон max=%.4f  баллистика max=%.4f рад/с\n', max(w_norm(1:nb)), max(w_norm(nb+6:end)));
fprintf('        порог скорости вращения Земли = %.3e рад/с\n', c.W_E);

% Сохраняем для следующего шага
save('falco_step12.mat', 'truth', 'imu', 'cfg', 'c');
fprintf('\nсохранено: falco_step12.mat\n');


% % === ВИЗУАЛИЗАЦИЯ ТРАЕКТОРИИ И ПРОФИЛЕЙ (вставить в конец MAIN_STEP12) ===
% load('falco_step12.mat');
% 
% dR = truth.R - truth.r0';
% alt = dR * truth.Cen(:,3);                    % высота (проекция на Up)
% rng = sqrt(dR(:,1).^2 + dR(:,2).^2);         % горизонтальная дальность
% 
% % === ФИГУРА 1: 3D ТРАЕКТОРИЯ В ECEF ===
% fig1 = figure('Color', 'w', 'Position', [100, 100, 950, 650]);
% 
% % Разделяем разгон и баллистику по цвету
% plot3(dR(1:nb,1)/1000, dR(1:nb,2)/1000, dR(1:nb,3)/1000, ...
%       'b-', 'LineWidth', 3, 'DisplayName', 'Разгон (0-10 с)');
% hold on;
% plot3(dR(nb:end,1)/1000, dR(nb:end,2)/1000, dR(nb:end,3)/1000, ...
%       'c--', 'LineWidth', 2.5, 'DisplayName', 'Баллистика');
% 
% % Ключевые точки
% plot3(0, 0, 0, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', ...
%       'DisplayName', 'Старт');
% [~, imax] = max(alt);
% plot3(dR(imax,1)/1000, dR(imax,2)/1000, dR(imax,3)/1000, ...
%       'ms', 'MarkerSize', 10, 'MarkerFaceColor', 'm', ...
%       'DisplayName', sprintf('Апогей (%.2f км)', alt(imax)/1000));
% plot3(dR(nb,1)/1000, dR(nb,2)/1000, dR(nb,3)/1000, ...
%       'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
%       'DisplayName', 'Выключение');
% plot3(dR(end,1)/1000, dR(end,2)/1000, dR(end,3)/1000, ...
%       'kx', 'MarkerSize', 12, 'LineWidth', 2, ...
%       'DisplayName', sprintf('Удар (t=%.1f с)', truth.t(end)));
% 
% grid on; box on;
% xlabel('\DeltaX в ECEF (км)', 'FontSize', 11);
% ylabel('\DeltaY в ECEF (км)', 'FontSize', 11);
% zlabel('\DeltaZ в ECEF (км)', 'FontSize', 11);
% title('Траектория Falco NAV в инерциальных осях ECEF', 'FontSize', 13, 'FontWeight', 'bold');
% legend('Location', 'best', 'FontSize', 10);
% axis equal vis3d;
% view(45, 30);
% rotate3d on;  % Интерактивное вращение
% 
% % === ФИГУРА 2: ПРОФИЛИ (высота, дальность, скорость, ускорение vs время) ===
% fig2 = figure('Color', 'w', 'Position', [1100, 100, 900, 700]);
% 
% % Подфигура 1: высота и дальность
% subplot(3,1,1);
% yyaxis left
% plot(truth.t, alt/1000, 'b-', 'LineWidth', 2.5);
% ylabel('Высота (км)', 'FontSize', 11, 'Color', 'b');
% ylim([min(alt)/1000, max(alt)/1000*1.1]);
% yyaxis right
% plot(truth.t, rng/1000, 'r-', 'LineWidth', 2.5);
% ylabel('Дальность (км)', 'FontSize', 11, 'Color', 'r');
% grid on; title('Геометрия траектории vs время', 'FontSize', 12, 'FontWeight', 'bold');
% xline(truth.t_b, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
% 
% % Подфигура 2: модуль скорости
% subplot(3,1,2);
% v_mag = vecnorm(truth.V, 2, 2);
% plot(truth.t, v_mag, 'g-', 'LineWidth', 2.5);
% ylabel('Скорость |V| (м/с)', 'FontSize', 11);
% grid on;
% xline(truth.t_b, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
% 
% % Подфигура 3: модуль кинематического ускорения
% subplot(3,1,3);
% a_mag = vecnorm(truth.A, 2, 2);
% semilogy(truth.t, a_mag, 'm-', 'LineWidth', 2.5);
% ylabel('Ускорение |A| (м/с²)', 'FontSize', 11);
% xlabel('Время (с)', 'FontSize', 11);
% grid on;
% xline(truth.t_b, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
% 
% % Общее оформление (все подфигуры)
% for i = 1:3
%     subplot(3,1,i);
%     set(gca, 'FontSize', 10);
% end
% 
% % === ФИГУРА 3: СКОРОСТНОЙ ПОРТРЕТ (ДАЛЬНОСТЬ vs ВЫСОТА, "дальномер") ===
% fig3 = figure('Color', 'w', 'Position', [100, 800, 700, 600]);
% plot(rng/1000, alt/1000, 'b-', 'LineWidth', 2.5);
% hold on;
% plot(0, 0, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
% plot(rng(imax)/1000, alt(imax)/1000, 'ms', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
% plot(rng(nb)/1000, alt(nb)/1000, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% plot(rng(end)/1000, alt(end)/1000, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
% grid on; box on;
% xlabel('Горизонтальная дальность (км)', 'FontSize', 11);
% ylabel('Высота (км)', 'FontSize', 11);
% title('Профиль траектории (Дальность × Высота)', 'FontSize', 13, 'FontWeight', 'bold');
% legend('Траектория', 'Старт', 'Апогей', 'Выключение', 'Удар', 'Location', 'best');
% axis equal;