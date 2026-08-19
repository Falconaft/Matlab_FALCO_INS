%% COMPARE_IMU_CANDIDATES  Шаг 4: сравнение кандидатов IMU по Монте-Карло
%
%  Прогоняет Монте-Карло чистой ИНС на 60-с коасте для всех профилей-кандидатов
%  и строит сравнительные графики. Требует falco_step12.mat (truth + imu_ideal).
%
%  ВНИМАНИЕ: это сравнение БЕЗ ESKF - показывает чистый инерциальный дрейф.
%  Финальное сравнение по КВО с оценкой bias фильтром - после Шага 5 (ESKF).
%
%  Настраиваемые параметры: N_mc (число прогонов), fnav, t_coast_start.

clear; clc;
load('falco_step12.mat', 'truth', 'imu', 'c');

% ===================== ПАРАМЕТРЫ =====================
N_mc          = 1;          % число прогонов Монте-Карло (МЕНЯЙ ПО ЖЕЛАНИЮ)
fnav          = 250;          % частота навигатора [Гц]
t_coast_start = 73.624;       % начало коаста [с]
base_seed     = 1000;         % базовый seed для воспроизводимости
cep_target    = 10;           % целевой КВО [м]

% ===================== СПИСОК КАНДИДАТОВ =====================
profiles = { imu_profile_stim300(), ...
             imu_profile_imup_s(),  ...
             imu_profile_pulse40(), ...
             imu_profile_pulse80() };
n_prof = numel(profiles);

% ===================== ПРОГОН МОНТЕ-КАРЛО =====================
results = cell(n_prof,1);
fprintf('=== СРАВНЕНИЕ КАНДИДАТОВ IMU (N_mc=%d, чистая ИНС, без ESKF) ===\n\n', N_mc);
fprintf('%-18s %10s %12s %10s %8s\n','Профиль','КВО[м]','медиана3D','R95[м]','<%d м',cep_target);
fprintf('%s\n', repmat('-',1,62));
for j = 1:n_prof
    mc = run_montecarlo_coast(profiles{j}, imu, truth, c, fnav, t_coast_start, N_mc, base_seed);
    results{j} = mc;
    frac_pass = 100*mean(mc.err_horiz < cep_target);
    fprintf('%-18s %9.1f %11.1f %9.1f %7.0f%%\n', ...
            mc.profile_name, mc.cep, median(mc.err_radial), mc.cep_r95, frac_pass);
end
fprintf('%s\n', repmat('-',1,62));

% ===================== ГРАФИКИ =====================
plot_imu_comparison(results, cep_target);

save('falco_step4_compare.mat', 'results', 'N_mc', 'fnav');
fprintf('\nсохранено: falco_step4_compare.mat\n');
