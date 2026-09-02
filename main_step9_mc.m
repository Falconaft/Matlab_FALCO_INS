%% MAIN_STEP9_MC  Шаг 9: Монте-Карло на РЕАЛИСТИЧНОЙ траектории
%
%  Новый базовый КВО для Pulse-40 на траектории с ненулевой удельной силой.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг 9.1: ПАРАМЕТРЫ КАМПАНИИ
% =====================================================================
N_mc       = 10;        % финальная кампания
base_seed  = 3001;
cep_target = 10;

% =====================================================================
% Шаг 9.2: КОНФИГУРАЦИЯ
% =====================================================================
% Базовая конфигурация; ниже — только специфика Монте-Карло.
cfg = falco_config();
% sigma_pos = 5 м, sigma_vel = 0.1 м/с, P0_pos = 5 м совпадают с базовыми.
% Частота записи диагностики задаётся ЯВНО в герцах, прореживание
% вычисляется из неё. Раньше стояло число 250, привязанное к f_nav=250
% (250 тактов = 1 с); при смене частоты навигатора смысл менялся бы молча.
cfg.f_diag     = 1;                                   % [Гц]
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);   % пересчёт после смены f_diag
cfg.cep_target = cep_target;

% Аномалия гравитационного поля (навигатор её не знает).
% Уклонение отвеса — ГЛАВНОЕ для КВО (горизонтальная ошибка):
% 5 угл.сек равнина, 15 холмистая местность, 30 горы.
cfg.defl_vert_sigma = 10 / 206265;      % [рад] СКО уклонения отвеса
cfg.grav_anom_sigma = 50e-5;            % [м/с²] СКО аномалии (50 мГал)

% --- B4: НЕКОМПЕНСИРОВАННОЕ РАССОГЛАСОВАНИЕ ЭПОХ GNSS/ИНС ---
% Измерение относится к t_meas = t_nav - offset, но фильтр обрабатывает
% его как текущее. Компенсация НЕ вводится — это часть базовой модели.
%
% Значение 2 мс выбрано как реалистичная верхняя граница для системы
% БЕЗ аппаратной синхронизации по PPS. По парному МК (test_b2b4_paired_mc)
% вклад составляет +0.75 м, то есть ~17% бюджета КВО; при 5 мс R95 уже
% выходит за целевые 10 м.
cfg.gnss_time_offset_enable = true;
cfg.gnss_time_offset        = 2e-3;     % [с]

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;
cfg.outage = [ 0.0            p.outage_boost_end;
               t_coast_start  t_end + 1.0 ];

% Выставка: СКО (разыгрывается на каждой реализации)
cfg.align_sigma = [0.5e-3; 0.5e-3; 4.0e-3];
cfg.P0_att      = 1.3 * cfg.align_sigma;
% =====================================================================
% НАЧАЛЬНАЯ ОШИБКА ПОЗИЦИИ И СКОРОСТИ (предстартовое усреднение GNSS)
% =====================================================================
% Изделие перед пуском неподвижно, решение GNSS усредняется. Остаточная
% погрешность и есть начальная ошибка навигатора. Прежде навигатор
% стартовал ТОЧНО из истины, что занижало ошибку на начальном участке.
%
% СОГЛАСОВАННОСТЬ: P0 задаётся ТЕМИ ЖЕ величинами, что и фактическая
% ошибка. Фильтр знает ровно ту неопределённость, которая у него есть —
% не переуверен и не недоуверен.
%
% ВАЖНО: это НЕ шум измерений GNSS в полёте. cfg.sigma_pos и
% cfg.sigma_vel остаются без изменений.
cfg.init_pos_sigma = 1.0*ones(3,1);     % [м]   1-sigma после усреднения
cfg.init_vel_sigma = 0.03*ones(3,1);    % [м/с] 1-sigma
cfg.P0_pos = cfg.init_pos_sigma;
cfg.P0_vel = cfg.init_vel_sigma;

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

fprintf('=== ШАГ 9: ФИНАЛЬНАЯ КАМПАНИЯ МОНТЕ-КАРЛО ===\n');
fprintf('%s\n', repmat('*',1,62));
fprintf('*  BASELINE: B2 (аномалия гравитации) + B4 (сдвиг GNSS %.1f мс)  *\n', ...
        cfg.gnss_time_offset*1e3);
fprintf('*  Компенсация сдвига НЕ вводится — это часть модели.           *\n');
fprintf('%s\n', repmat('*',1,62));
fprintf('Кандидат: %s\n', prof.name);
fprintf('Траектория: %.1f с, коаст %.1f с, окно GNSS %.1f с\n', ...
        t_end, p.coast_duration, t_coast_start - p.outage_boost_end);
fprintf('Выставка (СКО): курс %.1f мрад\n', cfg.align_sigma(3)*1e3);
fprintf('Реализаций: %d\n', N_mc);
fprintf('%s\n', repmat('-',1,62));
fprintf('СОСТАВ BASELINE:\n');
fprintf('  B2 аномалия поля : уклонение отвеса %.1f", аномалия %.0f мГал\n', ...
        cfg.defl_vert_sigma*206265, cfg.grav_anom_sigma*1e5);
fprintf('  B4 сдвиг GNSS    : %.1f мс (%s)\n', cfg.gnss_time_offset*1e3, ...
        ternary_str(cfg.gnss_time_offset_enable,'ВКЛЮЧЕН','выключен'));
fprintf('  плечо антенны    : [%.2f %.2f %.2f] м\n', cfg.lever);
fprintf('  шум GNSS в полёте: %.1f м / %.2f м/с при %.0f Гц\n', ...
        cfg.sigma_pos(1), cfg.sigma_vel(1), cfg.f_gnss);
fprintf('  нач. ошибка (СКО): %.2f м / %.3f м/с (предстартовое усреднение)\n', ...
        cfg.init_pos_sigma(1), cfg.init_vel_sigma(1));
fprintf('  калибровка IMU   : гиро %.3f, акс. %.3f\n', ...
        get_or(prof,'cal_factor_gyro',get_or(prof,'cal_factor',1)), ...
        get_or(prof,'cal_factor_accel',get_or(prof,'cal_factor',1)));
fprintf('%s\n\n', repmat('-',1,62));

% =====================================================================
% Шаг 9.3: ПРОГОН
% =====================================================================
% ВАЖНО: используется run_montecarlo_eskf2 (не v1). Только v2 вызывает
% draw_gravity_anomaly, поэтому при вызове v1 параметры cfg.defl_vert_sigma
% и cfg.grav_anom_sigma задавались, но НЕ ДЕЙСТВОВАЛИ — аномалия поля не
% применялась вовсе. Профиль передаётся дважды (истина и настройка фильтра
% совпадают); рассогласование профилей нужно только в бюджете ошибок.
mc = run_montecarlo_eskf2(prof, prof, imu, truth, c, cfg, N_mc, base_seed);

% =====================================================================
% Шаг 9.4: РЕЗУЛЬТАТ
% =====================================================================
fprintf('\n=== РЕЗУЛЬТАТ: %s ===\n', mc.prof_name);
fprintf('    BASELINE = B2 + B4(%.1f мс), N = %d\n', ...
        cfg.gnss_time_offset*1e3, N_mc);
fprintf('%s\n', repmat('=',1,54));
fprintf('  КВО50 (медиана)             : %8.2f м\n', mc.cep);
fprintf('  R95                         : %8.2f м\n', mc.r95);
fprintf('  Уложилось в %2.0f м            : %7.1f %%\n', cep_target, 100*mc.frac_pass);
fprintf('%s\n', repmat('-',1,54));
fprintf('  Систематика (ENU)           : [%+.2f %+.2f %+.2f] м\n', mc.bias_enu);
fprintf('%s\n', repmat('-',1,54));
% КАЧЕСТВО ОЦЕНКИ BIAS — пара (истинное полное смещение, остаток оценки)
% в момент перед коастом. Процент "сколько выучено" НЕ выводится: при
% малом истинном смещении отношение неустойчиво (наблюдались -536%).
g0  = 9.80665;
deg = pi/180;
fprintf('  BIAS ПЕРЕД КОАСТОМ (истинное полное / остаток оценки):\n');
fprintf('    акселерометр : %7.3f / %7.3f мг\n', ...
        mc.ba_true_pre_med/g0*1e3, mc.ba_resid_pre_med/g0*1e3);
fprintf('    гироскоп     : %7.3f / %7.3f °/ч\n', ...
        mc.bg_true_pre_med/deg*3600, mc.bg_resid_pre_med/deg*3600);
fprintf('%s\n', repmat('-',1,54));
fprintf('  ПЕРЕД ТЕРМИНАЛЬНЫМ КОАСТОМ (t = %.1f с):\n', mc.t_coast_start);
fprintf('    горизонтальная ошибка     : %8.3f м\n',    mc.pre_horiz_med);
fprintf('    |dv|                      : %8.4f м/с\n',  mc.pre_dv_med);
fprintf('    |dpsi|                    : %8.3f мрад\n', mc.pre_dpsi_med);
fprintf('%s\n', repmat('-',1,54));

% =====================================================================
% КАЧЕСТВО ОЦЕНКИ СМЕЩЕНИЯ ГИРОСКОПА, ПО ОСЯМ
% =====================================================================
% Основные показатели — RMS истинного смещения и RMSE остатка. Процент
% приведён справочно: при малом RMS он неустойчив.
% Осевая разбивка существенна: возбуждение по осям различается на порядки,
% и сводная норма скрывает, какая именно ось портит результат.
if isfield(mc,'gb')
    fprintf('  СМЕЩЕНИЕ ГИРОСКОПА (°/ч):\n');
    fprintf('    %-14s %8s %8s %8s %8s\n', '', 'X', 'Y', 'Z', 'норма');
    gb_line('RMS истина (пред)',  mc.gb.pre.rms_true_axis, mc.gb.pre.rms_true);
    gb_line('RMSE остаток(пред)', mc.gb.pre.rmse_axis,     mc.gb.pre.rmse);
    gb_line('RMS истина (кон.)',  mc.gb.end.rms_true_axis, mc.gb.end.rms_true);
    gb_line('RMSE остаток(кон.)', mc.gb.end.rmse_axis,     mc.gb.end.rmse);
    fprintf('    устранено, %%      : %7.1f %8.1f %8.1f %8.1f\n', ...
            mc.gb.end.eta_axis, mc.gb.end.eta);
    fprintf('    корр. b_g<->s_g   : %7.2f %8.2f %8.2f\n', mc.gb.rho_pre_med);
    fprintf('    (|корр.| -> 1 означает НЕРАЗДЕЛИМОСТЬ состояний:\n');
    fprintf('     dw = dbg + diag(w)*dsg, развести можно только\n');
    fprintf('     за счёт изменения w во времени)\n');
end
fprintf('%s\n', repmat('-',1,54));

% =====================================================================
% КАЧЕСТВО ОЦЕНКИ МАСШТАБНЫХ КОЭФФИЦИЕНТОВ
% =====================================================================
% Показатель устранённой ошибки:
%     eta = 100 * (1 - RMSE(s_est - s_true) / RMS(s_true))
% eta ~ 100%  оценка практически полная;
% eta ~   0%  оценка бесполезна (эквивалентна нулевой);
% eta <   0   оценка ВРЕДНА: внесено больше ошибки, чем убрано.
%
% Две точки: перед отключением GNSS и в конце. На коасте измерений нет,
% поэтому расхождение между точками показывает, сохраняется ли оценка.
if isfield(mc,'sf')
    fprintf('  МАСШТАБНЫЕ КОЭФФИЦИЕНТЫ (ppm):\n');
    fprintf('    %-14s %9s %9s %9s\n', '', 'истина', 'остаток', 'устранено');
    sf_line('гиро (перед)',  mc.sf.gyro_pre);
    sf_line('гиро (конец)',  mc.sf.gyro_end);
    sf_line('акс. (перед)',  mc.sf.accel_pre);
    sf_line('акс. (конец)',  mc.sf.accel_end);
    fprintf('    по осям XYZ, устранено %%:\n');
    fprintf('      гиро (конец) : %6.1f %6.1f %6.1f\n', mc.sf.gyro_end.eta_axis);
    fprintf('      акс. (конец) : %6.1f %6.1f %6.1f\n', mc.sf.accel_end.eta_axis);
end
fprintf('%s\n', repmat('-',1,54));

% =====================================================================
% СОСТОЯТЕЛЬНОСТЬ ФИЛЬТРА
% =====================================================================
% СТРОГИЙ ТЕСТ — это ANEES/ANIS: усреднение по АНСАМБЛЮ при фиксированном
% времени, где реализации независимы. Он на графике plot_consistency.
%
% Ниже приведены средние ПО ВРЕМЕНИ внутри реализаций. Они удобны как
% компактный индикатор уровня, но СТРОГИМ 95% ТЕСТОМ НЕ ЯВЛЯЮТСЯ:
% отсчёты внутри реализации коррелированы, поэтому chi-квадрат коридор
% к ним неприменим.
if isfield(mc,'nees_gnss_mean')
    fprintf('  СОСТОЯТЕЛЬНОСТЬ, средние по времени (идеал = %d):\n', mc.consist_dim);
    fprintf('    NEES6 в окне GNSS         : %8.2f\n', mc.nees_gnss_mean);
    fprintf('    NEES6 на коасте           : %8.2f\n', mc.nees_coast_mean);
    fprintf('    NIS  в окне GNSS          : %8.2f\n', mc.nis_tavg_mean);
    fprintf('    (строгий тест — ANEES/ANIS на графике состоятельности)\n');
end
if isfield(mc,'anees')
    frac_hi = 100*mean(mc.anees > mc.anees_hi, 'omitnan');
    fprintf('    ANEES выше коридора 95%%   : %7.1f %% времени\n', frac_hi);
end
if isfield(mc,'anis')
    frac_hi = 100*mean(mc.anis > mc.anis_hi, 'omitnan');
    fprintf('    ANIS  выше коридора 95%%   : %7.1f %% обновлений\n', frac_hi);
end
fprintf('%s\n', repmat('=',1,54));

% --- Фактически разыгранные аномалии гравитационного поля (B2) ---
if isfield(mc,'dg_info') && ~isempty(mc.dg_info)
    defl_tot = [mc.dg_info.defl_total_arcsec]';
    anom_v   = [mc.dg_info.anom_vert_mgal]';
    dg_h     = [mc.dg_info.dg_horiz_mps2]';
    fprintf('\n  РАЗЫГРАННАЯ АНОМАЛИЯ ПОЛЯ (задано %.1f", %.0f мГал):\n', ...
            cfg.defl_vert_sigma*206265, cfg.grav_anom_sigma*1e5);
    fprintf('    уклонение отвеса : медиана %.2f", СКО %.2f", макс %.2f"\n', ...
            median(defl_tot), std(defl_tot), max(defl_tot));
    fprintf('    аномалия вертик. : среднее %+.2f мГал, СКО %.2f мГал\n', ...
            mean(anom_v), std(anom_v));
    fprintf('    |dg| гориз.      : медиана %.3e м/с²\n', median(dg_h));
    fprintf('    ожидаемый вклад  : %.2f м за коаст %.0f с\n', ...
            0.5*median(dg_h)*p.coast_duration^2, p.coast_duration);
end

plot_montecarlo_cep(mc, cep_target);
plot_consistency(mc, p, cfg);
save('falco_step9_mc.mat', 'mc', 'cfg', 'prof', 'N_mc', 'base_seed');
fprintf('\nсохранено: falco_step9_mc.mat\n');


function s = ternary_str(cond, a, b)
%TERNARY_STR  Выбор строки по условию.
    if cond, s = a; else, s = b; end
end


function v = get_or(s, name, default)
%GET_OR  Чтение поля структуры со значением по умолчанию.
    if isfield(s, name), v = s.(name); else, v = default; end
end


function sf_line(nm, m)
%SF_LINE  Строка сводки по масштабным коэффициентам.
    fprintf('    %-14s %9.1f %9.1f %8.1f%%\n', nm, m.rms_true, m.rmse, m.eta);
end


function gb_line(nm, v_axis, v_norm)
%GB_LINE  Строка сводки по смещению гироскопа.
    fprintf('    %-14s %8.3f %8.3f %8.3f %8.3f\n', nm, v_axis, v_norm);
end
