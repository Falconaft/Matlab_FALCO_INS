%% TEST_B2B4_PAIRED_MC  Парный Монте-Карло: чувствительность к сдвигу GNSS
%
%  ЗАДАЧА: оценить, как некомпенсированное рассогласование эпох GNSS/ИНС
%  влияет на КВО в РЕАЛИСТИЧНЫХ условиях (полная модель ошибок Pulse-40,
%  включённая аномалия гравитации B2).
%
%  ПАРНОСТЬ — ключевое свойство теста. Для каждой реализации i все прогоны
%  при разных сдвигах используют ОДИН И ТОТ ЖЕ набор случайных величин:
%      seed_i           -> ошибки IMU и выставка
%      seed_i + OFFSET  -> аномалия гравитации (отдельный поток)
%  Различается ТОЛЬКО cfg.gnss_time_offset. Поэтому разность результатов
%  между сдвигами есть чистый эффект задержки, без вклада разброса выборки.
%  Это позволяет увидеть эффект на N = 25, где непарное сравнение
%  потребовало бы сотен реализаций.
%
%  ЧТО НЕ ДЕЛАЕТСЯ: компенсация задержки не вводится, архитектура фильтра
%  не меняется. Задача только измерить чувствительность.
%
%  Требует falco_step7.mat.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг P.1: ПАРАМЕТРЫ КАМПАНИИ
% =====================================================================
N_mc        = 25;
base_seed   = 9000;
offsets_ms  = [0, 1, 2, 5];
cep_target  = 10;

% Смещение потока для аномалии — то же, что в run_montecarlo_eskf2,
% чтобы реализации поля совпадали с обычными кампаниями Шага 9.
GRAV_SEED_OFFSET = 100000;

% =====================================================================
% Шаг P.2: КОНФИГУРАЦИЯ (как в актуальном Шаге 9)
% =====================================================================
cfg = falco_config();
cfg.f_diag     = 1;
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);
cfg.cep_target = cep_target;

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;
cfg.outage = [ 0.0            p.outage_boost_end;
               t_coast_start  t_end + 1.0 ];

cfg.align_sigma = [0.5e-3; 0.5e-3; 4.0e-3];
cfg.P0_att      = 1.3 * cfg.align_sigma;

% B2 ВКЛЮЧЕНА
cfg.defl_vert_sigma = 10 / 206265;      % уклонение отвеса, СКО [рад]
cfg.grav_anom_sigma = 50e-5;            % аномалия вертикали, СКО [м/с²]

prof = imu_profile_pulse40_updated();

fprintf('=== ПАРНЫЙ МК: B2 + B4 (чувствительность к сдвигу GNSS) ===\n');
fprintf('Кандидат: %s\n', prof.name);
fprintf('Реализаций: %d, базовый seed: %d\n', N_mc, base_seed);
fprintf('Сдвиги GNSS: %s мс\n', mat2str(offsets_ms));
fprintf('B2: уклонение отвеса %.1f", аномалия %.0f мГал\n', ...
        cfg.defl_vert_sigma*206265, cfg.grav_anom_sigma*1e5);
fprintf('Выставка (СКО): курс %.1f мрад\n', cfg.align_sigma(3)*1e3);
fprintf('Окно GNSS %.1f с, коаст %.1f с\n\n', ...
        t_coast_start - p.outage_boost_end, p.coast_duration);

% =====================================================================
% Шаг P.3: ПРОГОН — внешний цикл по сдвигам, внутренний по реализациям
% =====================================================================
n_off = numel(offsets_ms);

% Матрицы результатов: строка — реализация, столбец — сдвиг
E_horiz = nan(N_mc, n_off);   E_3d    = nan(N_mc, n_off);
E_dv    = nan(N_mc, n_off);
P_horiz = nan(N_mc, n_off);   P_dpsi  = nan(N_mc, n_off);
P_ba    = nan(N_mc, n_off);   P_dv    = nan(N_mc, n_off);
n_fail  = zeros(1, n_off);

t_start = tic;
for j = 1:n_off
    ms = offsets_ms(j);
    cfg_j = cfg;
    cfg_j.gnss_time_offset_enable = (ms > 0);
    cfg_j.gnss_time_offset        = ms * 1e-3;

    fprintf('--- сдвиг %.1f мс ---\n', ms);

    for i = 1:N_mc
        seed_i = base_seed + i;

        % ПАРНОСТЬ: те же потоки, что и при других сдвигах
        rng_imu  = RandStream('mt19937ar','Seed', seed_i);
        e_imu    = imu_draw_errors(prof, rng_imu);
        imu_meas = add_imu_errors(imu, prof, e_imu, rng_imu);

        cfg_i = cfg_j;
        cfg_i.align_err = cfg.align_sigma(:) .* randn(rng_imu, 3, 1);
        cfg_i.seed      = seed_i;

        rng_grav = RandStream('mt19937ar','Seed', seed_i + GRAV_SEED_OFFSET);
        cfg_i.dg_model = draw_gravity_anomaly(cfg, truth.Cen, rng_grav);

        res = eskf_run(imu_meas, truth, c, prof, cfg_i);

        if ~isfinite(res.err_horiz_final) || res.err_horiz_final > 1e5
            n_fail(j) = n_fail(j) + 1;
            continue;
        end

        E_horiz(i,j) = res.err_horiz_final;
        E_3d(i,j)    = res.err_3d_final;
        E_dv(i,j)    = norm(res.dv(end,:));

        % Состояние перед терминальным коастом
        [~, k_pre]   = min(abs(res.t - t_coast_start));
        P_horiz(i,j) = hypot(res.dr_ned(k_pre,1), res.dr_ned(k_pre,2));
        P_dv(i,j)    = norm(res.dv(k_pre,:));
        P_dpsi(i,j)  = norm(res.dpsi(k_pre,:)) * 1e3;
        P_ba(i,j)    = norm(res.ba_est(k_pre,:));

        if mod(i, max(1,round(N_mc/5))) == 0
            fprintf('    %3d/%3d, прошло %.0f с\n', i, N_mc, toc(t_start));
        end
    end
end
fprintf('\nКампания завершена за %.1f с (%.2f с/прогон)\n\n', ...
        toc(t_start), toc(t_start)/(N_mc*n_off));

% =====================================================================
% Шаг P.4: СВОДНАЯ ТАБЛИЦА
% =====================================================================
g0 = 9.80665;
fprintf('%s\n', repmat('=',1,104));
fprintf('  ПАРНЫЙ МК: N=%d на каждый сдвиг, B2 включена\n', N_mc);
fprintf('%s\n', repmat('=',1,104));
fprintf('  %7s | %8s %9s %8s %8s | %9s %10s | %9s %9s %9s %9s\n', ...
        'dt, мс', 'КВО50', 'среднее', 'R95', 'макс', '3D кон.', '|dv| кон.', ...
        'гориз.пред', '|dv|пред', 'dpsi пред', '|ba|пред');
fprintf('  %s\n', repmat('-',1,100));
for j = 1:n_off
    ok = ~isnan(E_horiz(:,j));
    fprintf('  %7.1f | %8.3f %9.3f %8.3f %8.3f | %9.3f %10.5f | %9.3f %9.5f %9.4f %9.5f\n', ...
        offsets_ms(j), median(E_horiz(ok,j)), mean(E_horiz(ok,j)), ...
        prctile(E_horiz(ok,j),95), max(E_horiz(ok,j)), median(E_3d(ok,j)), ...
        median(E_dv(ok,j)), ...
        median(P_horiz(ok,j)), median(P_dv(ok,j)), median(P_dpsi(ok,j)), ...
        median(P_ba(ok,j))/g0*1e3);
end
fprintf('%s\n', repmat('=',1,104));
fprintf('  "пред" — состояние непосредственно перед терминальным коастом (t=%.1f с)\n', ...
        t_coast_start);
fprintf('  |ba| в мг, dpsi в мрад, terminal velocity error в м/с\n');
fprintf('  Расходившихся прогонов: %s\n', mat2str(n_fail));

% =====================================================================
% Шаг P.5: ПАРНАЯ РАЗНОСТЬ (главный результат)
% =====================================================================
% Парность позволяет вычитать результаты ОДНОЙ И ТОЙ ЖЕ реализации,
% что убирает разброс выборки и оставляет чистый эффект сдвига.
fprintf('\n%s\n', repmat('=',1,84));
fprintf('  ПАРНАЯ РАЗНОСТЬ относительно dt = 0 (по каждой реализации)\n');
fprintf('%s\n', repmat('=',1,84));
fprintf('  %7s | %12s %12s %12s | %14s %12s\n', ...
        'dt, мс', 'd(гориз), м', 'СКО, м', 'на 1 мс', 'd(пред.гор), м', 'доля >0');
fprintf('  %s\n', repmat('-',1,80));
for j = 2:n_off
    d  = E_horiz(:,j) - E_horiz(:,1);
    dp = P_horiz(:,j) - P_horiz(:,1);
    ok = ~isnan(d);
    fprintf('  %7.1f | %12.4f %12.4f %12.4f | %14.4f %11.0f%%\n', ...
            offsets_ms(j), mean(d(ok)), std(d(ok)), ...
            mean(d(ok))/offsets_ms(j), mean(dp(ok)), 100*mean(d(ok)>0));
end
fprintf('%s\n', repmat('=',1,84));
fprintf('  "доля >0" — в скольких реализациях сдвиг УХУДШИЛ результат\n');

% --- Значимость парной разности ---
fprintf('\n  ЗНАЧИМОСТЬ (парный t-критерий против нулевой разности):\n');
for j = 2:n_off
    d  = E_horiz(:,j) - E_horiz(:,1);
    ok = ~isnan(d);
    se = std(d(ok))/sqrt(sum(ok));
    tt = mean(d(ok))/max(se,eps);
    fprintf('    dt=%4.1f мс: разность %+.4f м, СКО среднего %.4f, t = %+6.1f  [%s]\n', ...
            offsets_ms(j), mean(d(ok)), se, tt, ...
            ternary_str(abs(tt) > 2, 'ЗНАЧИМО', 'в пределах шума'));
end

% =====================================================================
% Шаг P.6: ГРАФИКИ
% =====================================================================
off = offsets_ms(:);
figure('Color','w','Position',[60 60 1100 660]);

ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
cep = arrayfun(@(j) median(E_horiz(~isnan(E_horiz(:,j)),j)), 1:n_off)';
r95 = arrayfun(@(j) prctile(E_horiz(~isnan(E_horiz(:,j)),j),95), 1:n_off)';
mx  = arrayfun(@(j) max(E_horiz(~isnan(E_horiz(:,j)),j)), 1:n_off)';
plot(ax, off, cep, '-o', 'LineWidth',1.9);
plot(ax, off, r95, '-s', 'LineWidth',1.5);
plot(ax, off, mx,  '-^', 'LineWidth',1.2);
yline(ax, cep_target, 'k--', 'LineWidth',1.4);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'Горизонтальная ошибка, м');
title(ax,'КВО и хвосты распределения','FontWeight','bold');
legend(ax,{'КВО50','R95','макс','цель'},'Location','best');

ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for j = 1:n_off
    d = E_horiz(:,j) - E_horiz(:,1);
    plot(ax, offsets_ms(j)*ones(N_mc,1), d, 'o', 'MarkerSize',4);
end
dm = arrayfun(@(j) mean(E_horiz(~isnan(E_horiz(:,j)),j) - E_horiz(~isnan(E_horiz(:,j)),1)), 1:n_off)';
plot(ax, off, dm, '-k', 'LineWidth',2);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'Парная разность, м');
title(ax,'Парная разность по реализациям','FontWeight','bold');

ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
pre = arrayfun(@(j) median(P_horiz(~isnan(P_horiz(:,j)),j)), 1:n_off)';
plot(ax, off, pre, '-o', 'LineWidth',1.8);
plot(ax, off, cep, '-s', 'LineWidth',1.5);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'Горизонтальная ошибка, м');
title(ax,'Перед коастом и на ударе','FontWeight','bold');
legend(ax,{'перед коастом','на ударе'},'Location','best');

ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
ba = arrayfun(@(j) median(P_ba(~isnan(P_ba(:,j)),j)), 1:n_off)'/g0*1e3;
plot(ax, off, ba, '-o', 'LineWidth',1.8, 'Color',[0.8 0.3 0.2]);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'|ba| перед коастом, мг');
title(ax,'Поглощается ли сдвиг в оценку bias','FontWeight','bold');

save('falco_test_b2b4_mc.mat', 'E_horiz','E_3d','E_dv', ...
     'P_horiz','P_dpsi','P_ba','P_dv','offsets_ms','N_mc','base_seed','cfg');
fprintf('\nсохранено: falco_test_b2b4_mc.mat\n');


function s = ternary_str(cond, a, b)
%TERNARY_STR  Выбор строки по условию.
    if cond, s = a; else, s = b; end
end
