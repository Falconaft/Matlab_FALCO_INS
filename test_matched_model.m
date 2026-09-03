%% TEST_MATCHED_MODEL  Сравнение baseline и согласованной модели
%
%  ЗАДАЧА. Baseline показал: фильтр СОСТОЯТЕЛЕН по NIS, но НЕ состоятелен
%  по NEES. Гипотеза — причина в источниках, которые ЕСТЬ в истине, но
%  ОТСУТСТВУЮТ в модели фильтра (нет ни состояний, ни вклада в Q).
%
%  ПРОВЕРКА. Прогоняются две кампании на ОДНИХ И ТЕХ ЖЕ seeds:
%
%    baseline        — истина содержит всё, фильтр как есть;
%    matched-model   — из ИСТИНЫ убраны источники, не представленные
%                      в модели фильтра; ФИЛЬТР НЕ МЕНЯЕТСЯ ВОВСЕ.
%
%  Отключаются в истине:
%    B2 аномалия гравитации, B4 сдвиг GNSS, перекос осей обеих триад,
%    остаточная g-чувствительность гироскопа, VRC обоих датчиков.
%
%  Остаются в истине (они ПРЕДСТАВЛЕНЫ в состоянии фильтра):
%    turn-on и in-run bias, scale factor, белый шум ARW/VRW, выставка.
%
%  ВАЖНО: P, Q, R, F, H и состав состояния НЕ трогаются. prof_filter
%  передаётся исходным. Используется штатный run_montecarlo_eskf2
%  с раздельными профилями истины и фильтра.
%
%  Требует falco_step7.mat.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг M.1: ПАРАМЕТРЫ
% =====================================================================
N_mc      = 20;
base_seed = 3001;

% =====================================================================
% Шаг M.2: КОНФИГУРАЦИЯ BASELINE (как в main_step9_mc)
% =====================================================================
cfg = falco_config();
cfg.f_diag     = 10;
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);
cfg.cep_target = 10;

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;
cfg.outage = [ 0.0            p.outage_boost_end;
               t_coast_start  t_end + 1.0 ];

cfg.align_sigma = [0.5e-3; 0.5e-3; 4.0e-3];
cfg.P0_att      = 1.3 * cfg.align_sigma;

cfg.init_pos_sigma = 1.0*ones(3,1);
cfg.init_vel_sigma = 0.03*ones(3,1);
cfg.P0_pos = cfg.init_pos_sigma;
cfg.P0_vel = cfg.init_vel_sigma;

cfg.defl_vert_sigma = 10 / 206265;
cfg.grav_anom_sigma = 50e-5;
cfg.gnss_time_offset_enable = true;
cfg.gnss_time_offset        = 2e-3;

prof = imu_profile_pulse40_updated();

fprintf('=== MATCHED-MODEL TEST ===\n');
fprintf('N = %d, base_seed = %d, кандидат %s\n', N_mc, base_seed, prof.name);
fprintf('Фильтр НЕ меняется ни в одной кампании (P/Q/R/F/H исходные).\n\n');

% =====================================================================
% Шаг M.3: КАМПАНИЯ 1 — BASELINE
% =====================================================================
fprintf('--- КАМПАНИЯ 1: baseline (истина содержит всё) ---\n');
mc_base = run_montecarlo_eskf2(prof, prof, imu, truth, c, cfg, N_mc, base_seed);

% =====================================================================
% Шаг M.4: КАМПАНИЯ 2 — MATCHED MODEL
% =====================================================================
% Из ИСТИНЫ убираются источники, не представленные в модели фильтра.
% Копии структур; рабочие prof и cfg не трогаются.
prof_m = prof;
prof_m.accel.misalign_axes  = 0*prof_m.accel.misalign_axes;
prof_m.accel.misalign_ortho = 0*prof_m.accel.misalign_ortho;
prof_m.gyro.misalign_axes   = 0*prof_m.gyro.misalign_axes;
prof_m.gyro.misalign_ortho  = 0*prof_m.gyro.misalign_ortho;

prof_m.gyro.g_sens_cal_sigma = 0;   % остаточная g-чувствительность
prof_m.accel.vrc = 0;               prof_m.gyro.vrc = 0;

cfg_m = cfg;
cfg_m.defl_vert_sigma = 0;          % B2 выкл
cfg_m.grav_anom_sigma = 0;
cfg_m.gnss_time_offset_enable = false;   % B4 выкл
cfg_m.gnss_time_offset        = 0;

fprintf('\n--- КАМПАНИЯ 2: matched-model ---\n');
fprintf('Отключено в ИСТИНЕ: B2, B4, перекос осей (акс+гиро),\n');
fprintf('остаток g-sens, VRC (акс+гиро). Фильтр — prof исходный.\n');
mc_match = run_montecarlo_eskf2(prof_m, prof, imu, truth, c, cfg_m, N_mc, base_seed);

% =====================================================================
% Шаг M.4b: КОНТРОЛЬ ПАРНОСТИ РЕАЛИЗАЦИЙ
% =====================================================================
% Обе кампании идут по одному base_seed. Розыгрыш выставки и начальной
% ошибки происходит из потока rng_imu ПОСЛЕ ошибок IMU, а обнуление
% параметров профиля НЕ меняет числа обращений к randn (поля остаются,
% меняются лишь значения). Поэтому эти величины ОБЯЗАНЫ совпасть
% побитово. Расхождение означало бы рассинхронизацию потоков, при которой
% сравнение кампаний теряет смысл.
fprintf('\n--- КОНТРОЛЬ ПАРНОСТИ ---\n');
chk = { 'align_used',    'выставка';
        'init_pos_used', 'нач. ошибка позиции';
        'init_vel_used', 'нач. ошибка скорости' };
all_ok = true;
for q = 1:size(chk,1)
    f = chk{q,1};
    if isfield(mc_base,f) && isfield(mc_match,f)
        d = max(abs(mc_base.(f)(:) - mc_match.(f)(:)));
        ok = (d == 0);
        all_ok = all_ok && ok;
        fprintf('  %-24s макс |разность| = %.3e  [%s]\n', chk{q,2}, d, ...
                pass_str(ok));
    else
        fprintf('  %-24s поле отсутствует\n', chk{q,2});
        all_ok = false;
    end
end
if all_ok
    fprintf('  Реализации попарно идентичны — сравнение корректно.\n');
else
    warning('test_matched_model:pairingBroken', ...
        ['Реализации НЕ идентичны. Потоки случайных чисел разошлись, ' ...
         'сравнение кампаний некорректно.']);
end

% =====================================================================
% Шаг M.5: СРАВНЕНИЕ
% =====================================================================
fprintf('\n%s\n', repmat('=',1,72));
fprintf('  СРАВНЕНИЕ: baseline против matched-model (N = %d)\n', N_mc);
fprintf('%s\n', repmat('=',1,72));
fprintf('  %-28s %14s %14s %10s\n', 'Метрика', 'baseline', 'matched', 'отношение');
fprintf('  %s\n', repmat('-',1,68));

pr = @(nm, a, b) fprintf('  %-28s %14.3f %14.3f %10.2f\n', nm, a, b, a/max(b,eps));

pr('КВО50, м',            mc_base.cep,  mc_match.cep);
pr('R95, м',              mc_base.r95,  mc_match.r95);

fprintf('  %s\n', repmat('-',1,68));
fprintf('  СОСТОЯТЕЛЬНОСТЬ (средние ANEES/ANIS по времени; идеал в скобках)\n');
mB = @(f) mean(mc_base.(f)(isfinite(mc_base.(f))));
mM = @(f) mean(mc_match.(f)(isfinite(mc_match.(f))));
pr('ANEES6      (идеал 6)',  mB('anees'),   mM('anees'));
if isfield(mc_base,'anees_r')
    pr('ANEES_r     (идеал 3)', mB('anees_r'), mM('anees_r'));
    pr('ANEES_v     (идеал 3)', mB('anees_v'), mM('anees_v'));
end
pr('ANIS        (идеал 6)',  mB('anis'),    mM('anis'));

fprintf('  %s\n', repmat('-',1,68));
fprintf('  ДОЛЯ ВРЕМЕНИ ВЫШЕ КОРИДОРА 95%%\n');
fo = @(m,f,fh) 100*mean(m.(f) > m.(fh), 'omitnan');
fprintf('  %-28s %13.1f%% %13.1f%%\n', 'ANEES6', ...
        fo(mc_base,'anees','anees_hi'), fo(mc_match,'anees','anees_hi'));
if isfield(mc_base,'anees_r')
    fprintf('  %-28s %13.1f%% %13.1f%%\n', 'ANEES_r', ...
            fo(mc_base,'anees_r','anees_r_hi'), fo(mc_match,'anees_r','anees_r_hi'));
    fprintf('  %-28s %13.1f%% %13.1f%%\n', 'ANEES_v', ...
            fo(mc_base,'anees_v','anees_v_hi'), fo(mc_match,'anees_v','anees_v_hi'));
end
fprintf('  %-28s %13.1f%% %13.1f%%\n', 'ANIS', ...
        fo(mc_base,'anis','anis_hi'), fo(mc_match,'anis','anis_hi'));

if isfield(mc_base,'rho_pos_med')
    fprintf('  %s\n', repmat('-',1,68));
    fprintf('  ВКЛАД P ОТНОСИТЕЛЬНО R: rho = diag(H·P_prior·H'')/diag(R)\n');
    fprintf('  (извлечено из той же S, что использует NIS)\n');
    pr('rho позиция (медиана)', mc_base.rho_pos_med, mc_match.rho_pos_med);
    pr('rho скорость (медиана)', mc_base.rho_vel_med, mc_match.rho_vel_med);
    fprintf('  При rho << 1 величина S определяется шумом R, и NIS\n');
    fprintf('  теряет чувствительность к заниженности P.\n');
end
fprintf('%s\n', repmat('=',1,72));

% =====================================================================
% Шаг M.6: ГРАФИКИ
% =====================================================================
plot_matched_comparison(mc_base, mc_match, cfg);

save('falco_test_matched.mat', 'mc_base', 'mc_match', 'cfg', 'cfg_m', ...
     'prof', 'prof_m', 'N_mc', 'base_seed');
fprintf('\nсохранено: falco_test_matched.mat\n');


function s = pass_str(ok)
%PASS_STR  Метка результата проверки.
    if ok, s = 'OK'; else, s = 'РАСХОЖДЕНИЕ'; end
end
