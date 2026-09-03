%% TEST_Q_MODEL_SWEEP  Короткий sweep модельного шума Q на realistic baseline
%
%  ЗАДАЧА. Matched-model тест показал, что ядро ESKF статистически
%  состоятельно, но realistic baseline имеет сильно завышенный NEES из-за
%  известных источников model mismatch.
%
%  Здесь проверяется простой инженерный tuning:
%  умеренно увеличить process noise фильтра, НЕ меняя truth, P0, R, F, H
%  и размерность состояния.
%
%  Для минимального вмешательства используется существующая возможность
%  run_montecarlo_eskf2 передавать РАЗНЫЕ профили истины и фильтра:
%
%    prof_truth  = исходный Pulse-40 во всех кампаниях;
%    prof_filter = копия Pulse-40, где ARW/VRW умножены на sqrt(kQ).
%
%  Поскольку eskf_Q_matrix использует ARW^2 и VRW^2:
%
%       Q_white(filter) = kQ * Q_white(original)
%
%  При этом GM process noise bias и scale-factor Q НЕ меняются.
%
%  Требует falco_step7.mat.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг Q.1: ПАРАМЕТРЫ SWEEP
% =====================================================================
N_mc      = 10;
base_seed = 3001;

% Множитель только white process Q:
% 1 = штатный Q, 2 = удвоенный, 4 = учетверённый
kQ_list = [1 2 4];

% =====================================================================
% Шаг Q.2: REALISTIC BASELINE — как в main_step9_mc
% =====================================================================
cfg = falco_config();

cfg.f_diag     = 10;  % детальная диагностика для tuning
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);
cfg.cep_target = 10;

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;

cfg.outage = [ 0.0            p.outage_boost_end;
               t_coast_start  t_end + 1.0 ];

% Выставка
cfg.align_sigma = [0.5e-3; 0.5e-3; 4.0e-3];
cfg.P0_att      = 1.3 * cfg.align_sigma;

% Начальные ошибки после предстартового усреднения GNSS
cfg.init_pos_sigma = 1.0*ones(3,1);
cfg.init_vel_sigma = 0.03*ones(3,1);

cfg.P0_pos = cfg.init_pos_sigma;
cfg.P0_vel = cfg.init_vel_sigma;

% B2 — аномалия гравитационного поля
cfg.defl_vert_sigma = 10 / 206265;
cfg.grav_anom_sigma = 50e-5;

% B4 — некомпенсированное рассогласование эпох GNSS/INS
cfg.gnss_time_offset_enable = true;
cfg.gnss_time_offset        = 2e-3;

prof_truth = imu_profile_pulse40_updated();

fprintf('=== Q-MODEL SWEEP TEST ===\n');
fprintf('N = %d на вариант, base_seed = %d, кандидат %s\n', ...
        N_mc, base_seed, prof_truth.name);
fprintf('Truth неизменна во всех кампаниях.\n');
fprintf('P0/R/F/H неизменны.\n');
fprintf('Меняется только white process Q фильтра через ARW/VRW.\n');
fprintf('kQ = %s\n\n', mat2str(kQ_list));

% =====================================================================
% Шаг Q.3: КАМПАНИИ
% =====================================================================
mc = cell(numel(kQ_list),1);
prof_filter_all = cell(numel(kQ_list),1);

for iq = 1:numel(kQ_list)

    kQ = kQ_list(iq);

    % eskf_Q_matrix использует ARW^2 и VRW^2.
    % Поэтому sqrt(kQ) по коэффициенту -> kQ по мощности Q.
    prof_filter = prof_truth;

    prof_filter.gyro.arw = ...
        prof_truth.gyro.arw * sqrt(kQ);

    prof_filter.accel.vrw = ...
        prof_truth.accel.vrw * sqrt(kQ);

    prof_filter_all{iq} = prof_filter;

    fprintf('--- КАМПАНИЯ %d/%d: kQ = %.3g ---\n', ...
            iq, numel(kQ_list), kQ);

    mc{iq} = run_montecarlo_eskf2( ...
        prof_truth, ...
        prof_filter, ...
        imu, truth, c, cfg, ...
        N_mc, base_seed);

    fprintf('  CEP50 = %.3f м, R95 = %.3f м\n\n', ...
            mc{iq}.cep, mc{iq}.r95);
end

% =====================================================================
% Шаг Q.4: КОНТРОЛЬ ПАРНОСТИ
% =====================================================================
fprintf('--- КОНТРОЛЬ ПАРНОСТИ ---\n');

chk = { ...
    'align_used',    'выставка';
    'init_pos_used', 'нач. ошибка позиции';
    'init_vel_used', 'нач. ошибка скорости' };

all_ok = true;

for iq = 2:numel(kQ_list)

    fprintf('  kQ=%.3g против kQ=%.3g:\n', ...
            kQ_list(iq), kQ_list(1));

    for q = 1:size(chk,1)

        f = chk{q,1};

        d = max(abs(mc{1}.(f)(:) - mc{iq}.(f)(:)));

        ok = (d == 0);
        all_ok = all_ok && ok;

        fprintf('    %-22s max |diff| = %.3e  [%s]\n', ...
                chk{q,2}, d, pass_str(ok));
    end
end

if all_ok
    fprintf('  Все кампании попарно идентичны.\n');
else
    warning('test_q_model_sweep:pairingBroken', ...
        'Парность кампаний нарушена.');
end

% =====================================================================
% Шаг Q.5: СВОДНАЯ СТАТИСТИКА
% =====================================================================
K = numel(kQ_list);

cep       = nan(K,1);
r95       = nan(K,1);
pass10    = nan(K,1);

anees6    = nan(K,1);
anees_r   = nan(K,1);
anees_v   = nan(K,1);
anis      = nan(K,1);

frac_n6hi = nan(K,1);
frac_nrhi = nan(K,1);
frac_nvhi = nan(K,1);
frac_ishi = nan(K,1);

rho_pos   = nan(K,1);
rho_vel   = nan(K,1);

for iq = 1:K

    m = mc{iq};

    cep(iq)    = m.cep;
    r95(iq)    = m.r95;
    pass10(iq) = 100*m.frac_pass;

    anees6(iq) = finite_mean(m.anees);
    anees_r(iq)= finite_mean(m.anees_r);
    anees_v(iq)= finite_mean(m.anees_v);
    anis(iq)   = finite_mean(m.anis);

    frac_n6hi(iq) = ...
        100*mean(m.anees > m.anees_hi, 'omitnan');

    frac_nrhi(iq) = ...
        100*mean(m.anees_r > m.anees_r_hi, 'omitnan');

    frac_nvhi(iq) = ...
        100*mean(m.anees_v > m.anees_v_hi, 'omitnan');

    frac_ishi(iq) = ...
        100*mean(m.anis > m.anis_hi, 'omitnan');

    if isfield(m,'rho_pos_med')
        rho_pos(iq) = m.rho_pos_med;
    end

    if isfield(m,'rho_vel_med')
        rho_vel(iq) = m.rho_vel_med;
    end
end

% =====================================================================
% Шаг Q.6: ПЕЧАТЬ ТАБЛИЦЫ
% =====================================================================
fprintf('\n%s\n', repmat('=',1,100));
fprintf('Q-MODEL SWEEP: realistic baseline, N = %d\n', N_mc);
fprintf('%s\n', repmat('=',1,100));

fprintf('%-8s %9s %9s %9s %9s %9s %9s %9s\n', ...
    'kQ', 'CEP50', 'R95', '<10m,%', ...
    'ANEES6', 'ANEES_r', 'ANEES_v', 'ANIS');

fprintf('%s\n', repmat('-',1,100));

for iq = 1:K
    fprintf('%-8.3g %9.3f %9.3f %9.1f %9.3f %9.3f %9.3f %9.3f\n', ...
        kQ_list(iq), ...
        cep(iq), r95(iq), pass10(iq), ...
        anees6(iq), anees_r(iq), anees_v(iq), anis(iq));
end

fprintf('%s\n', repmat('-',1,100));

fprintf('ДОЛЯ ВРЕМЕНИ ВЫШЕ ВЕРХНЕЙ ГРАНИЦЫ 95%%\n');

fprintf('%-8s %12s %12s %12s %12s\n', ...
    'kQ', 'ANEES6,%', 'ANEES_r,%', 'ANEES_v,%', 'ANIS,%');

for iq = 1:K
    fprintf('%-8.3g %12.1f %12.1f %12.1f %12.1f\n', ...
        kQ_list(iq), ...
        frac_n6hi(iq), ...
        frac_nrhi(iq), ...
        frac_nvhi(iq), ...
        frac_ishi(iq));
end

fprintf('%s\n', repmat('-',1,100));

fprintf('rho = diag(H*P_prior*H'') / diag(R), медиана\n');

fprintf('%-8s %14s %14s\n', ...
    'kQ', 'rho position', 'rho velocity');

for iq = 1:K
    fprintf('%-8.3g %14.5f %14.5f\n', ...
        kQ_list(iq), ...
        rho_pos(iq), rho_vel(iq));
end

fprintf('%s\n', repmat('=',1,100));

% =====================================================================
% Шаг Q.7: ГРАФИК ANEES6
% =====================================================================
figure('Name','Q sweep - ANEES6','Color','w');

hold on;
grid on;

for iq = 1:K
    plot( ...
        mc{iq}.anees_t, ...
        mc{iq}.anees, ...
        'LineWidth',1.3, ...
        'DisplayName',sprintf('kQ = %.3g',kQ_list(iq)));
end

yline(6,'k-','идеал 6','LineWidth',1.2);
xline(t_coast_start,'k--','GNSS off');

set(gca,'YScale','log');

xlabel('Время, с');
ylabel('ANEES6');

title(sprintf('ANEES6(t): Q sweep, N=%d',N_mc));

legend('Location','best');

% =====================================================================
% Шаг Q.8: ГРАФИК ANIS
% =====================================================================
figure('Name','Q sweep - ANIS','Color','w');

hold on;
grid on;

for iq = 1:K
    plot( ...
        mc{iq}.anis_t, ...
        mc{iq}.anis, ...
        'LineWidth',1.0, ...
        'DisplayName',sprintf('kQ = %.3g',kQ_list(iq)));
end

yline(6,'k-','идеал 6','LineWidth',1.2);

xlabel('Время, с');
ylabel('ANIS');

title(sprintf('ANIS(t): Q sweep, N=%d',N_mc));

legend('Location','best');

% =====================================================================
% Шаг Q.9: ИНЖЕНЕРНАЯ СВОДКА
% =====================================================================
figure('Name','Q sweep - summary','Color','w');

tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

nexttile;

plot(kQ_list,anees6,'-o','LineWidth',1.3);
hold on;
grid on;

plot(kQ_list,anees_r,'-o','LineWidth',1.3);
plot(kQ_list,anees_v,'-o','LineWidth',1.3);

yline(6,'k--','ideal NEES6');
yline(3,'k:','ideal channel');

xlabel('kQ');
ylabel('mean ANEES');

title('Consistency vs Q');

legend( ...
    'ANEES6', ...
    'ANEES_r', ...
    'ANEES_v', ...
    'Location','best');

nexttile;

yyaxis left;

plot(kQ_list,cep,'-o','LineWidth',1.3);
hold on;
grid on;

plot(kQ_list,r95,'-o','LineWidth',1.3);

ylabel('Ошибка, м');

yyaxis right;

plot(kQ_list,pass10,'-o','LineWidth',1.3);

ylabel('<10 м, %');

xlabel('kQ');

title('Navigation quality vs Q');

legend( ...
    'CEP50', ...
    'R95', ...
    '<10 m', ...
    'Location','best');

% =====================================================================
% Шаг Q.10: СОХРАНЕНИЕ
% =====================================================================
summary = table( ...
    kQ_list(:), ...
    cep, r95, pass10, ...
    anees6, anees_r, anees_v, anis, ...
    frac_n6hi, frac_nrhi, frac_nvhi, frac_ishi, ...
    rho_pos, rho_vel, ...
    'VariableNames', { ...
        'kQ', ...
        'CEP50_m', ...
        'R95_m', ...
        'Pass10_pct', ...
        'ANEES6', ...
        'ANEES_r', ...
        'ANEES_v', ...
        'ANIS', ...
        'ANEES6_hi_pct', ...
        'ANEES_r_hi_pct', ...
        'ANEES_v_hi_pct', ...
        'ANIS_hi_pct', ...
        'rho_pos_med', ...
        'rho_vel_med' });

save('falco_test_q_model_sweep.mat', ...
    'mc', ...
    'summary', ...
    'cfg', ...
    'prof_truth', ...
    'prof_filter_all', ...
    'N_mc', ...
    'base_seed', ...
    'kQ_list');

fprintf('\nсохранено: falco_test_q_model_sweep.mat\n');


function y = finite_mean(x)
%FINITE_MEAN  Среднее только по конечным значениям.

    x = x(isfinite(x));

    if isempty(x)
        y = NaN;
    else
        y = mean(x);
    end
end


function s = pass_str(ok)
%PASS_STR  Метка результата проверки.

    if ok
        s = 'OK';
    else
        s = 'РАСХОЖДЕНИЕ';
    end
end