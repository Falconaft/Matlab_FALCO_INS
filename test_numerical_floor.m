%% TEST_NUMERICAL_FLOOR  Измерение собственного численного пола реализации
%
%  ДВА НЕЗАВИСИМЫХ ТЕСТА, которые НЕЛЬЗЯ смешивать:
%
%  TEST A — MECHANIZATION FLOOR
%     Идеальная ИНС, GNSS ОТКЛЮЧЁН НА ВСЁМ ПОЛЁТЕ.
%     Причина отключения: идеальный GNSS (sigma -> 0) корректировал бы
%     численную ошибку механизации и тем самым МАСКИРОВАЛ бы именно ту
%     величину, которую мы хотим измерить. Пол механизации измерим только
%     в разомкнутом контуре.
%     Измеряет: накопление ошибки чисто от дискретизации (RK4, трапеции
%     в приращениях, поворотная поправка, конинг/скаллинг, дискретизация
%     кватерниона), при нулевых ошибках датчиков и нулевых начальных ошибках.
%
%  TEST B — FULL LOOP IDEAL
%     Идеальная IMU + почти идеальный GNSS + работающий ESKF.
%     Измеряет пол ЗАМКНУТОГО контура: сколько остаётся после коррекции.
%     Ожидается заметно меньше, чем в TEST A, на участках с GNSS,
%     и рост на терминальном коасте.
%
%  ЧТО ОБНУЛЯЕТСЯ (только в КОПИЯХ структур, рабочие prof/cfg не трогаются):
%     turn-on bias, in-run bias, ARW, VRW, scale factor, misalignment,
%     orthogonality, ошибка калибровки g-sensitivity, VRC, вибрация,
%     начальные ошибки позиции/скорости/ориентации, аномалия гравитации.
%
%  ПРИМЕЧАНИЕ о g-sensitivity: обнуляется ОШИБКА КАЛИБРОВКИ, а сам
%  коэффициент g_sensitivity остаётся. Компенсация в механизации при этом
%  точна с точностью до дискретизации, и её остаток ВХОДИТ в измеряемый
%  пол — это честно, поскольку компенсация является частью механизации.
%
%  ОПОРНЫЕ ЗНАЧЕНИЯ BASELINE_250 (f_truth = 1000 Гц, f_nav = 250 Гц),
%  TEST A — mechanization floor:
%      terminal |dr|   = 1.512289e-02 м
%      terminal |dv|   = 1.303547e-04 м/с
%      terminal |dpsi| = 3.597268e-06 мрад
%      max |dr|        = 1.512289e-02 м (в t = 208.57 с)
%      эпох            = 52143
%  Критерий регрессии: тот же порядок величин или улучшение. Точного
%  совпадения не требуется — меняются и сетка истины, и частота навигатора.
%
%  Требует falco_step7.mat (результат main_step7), ПЕРЕГЕНЕРИРОВАННЫЙ после
%  изменения p.dt в traj_params.m. При старом файле (truth.dt = 0.001)
%  make_increments прервёт выполнение по assert, так как 1000/400 не целое.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг T.1: БАЗОВАЯ КОНФИГУРАЦИЯ (та же, что в рабочем сценарии Шага 8)
% =====================================================================
% ЧАСТОТА НАВИГАТОРА — главный параметр регрессионного сравнения.
%   400 Гц — режим Шага B (сетка истины 2000 Гц, 2000/400 = 5, кратно).
%   250 Гц — воспроизводит режим BASELINE_250; при сетке 2000 Гц отношение
%            2000/250 = 8 тоже целое, поэтому этим значением можно отделить
%            влияние более мелкой сетки истины от влияния самой частоты.
% Базовая конфигурация; ниже — специфика теста.
cfg = falco_config();
cfg.seed   = 3000;
cfg.f_diag = 250;                  % частота записи диагностики [Гц]
% ПРИМЕЧАНИЕ: diag_decim здесь НЕ пересчитывается специально — ниже, при
% формировании идеальных копий, он принудительно ставится в 1, чтобы
% диагностика писалась на каждом такте и максимумы не терялись.

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;

cfg.outage = [ 0.0            p.outage_boost_end;
               t_coast_start  t_end + 1.0 ];

cfg.align_err = [0.5e-3; 0.5e-3; 4e-3];
cfg.P0_att    = 1.3 * abs(cfg.align_err);
cfg.P0_pos    = [5.0; 5.0; 5.0];
cfg.P0_vel    = [0.1; 0.1; 0.1];

prof = imu_profile_pulse40_updated();

% =====================================================================
% Шаг T.2: ПРОВЕРКА ВРЕМЕННОЙ СЕТКИ (A1)
% =====================================================================
print_time_grid(truth, cfg.fnav, cfg.f_gnss);

% =====================================================================
% Шаг T.2b: РЕГРЕССИОННЫЕ ПРОВЕРКИ NED-ИНТЕРФЕЙСА
% =====================================================================
% Проверки вынесены сюда намеренно: C_e_ned вызывается в циклах, поэтому
% контроль ортогональности внутри неё был бы лишней работой на каждом шаге.
run_ned_checks(truth);

% =====================================================================
% Шаг T.3: ПОДГОТОВКА ИДЕАЛЬНЫХ КОПИЙ
% =====================================================================
% ВАЖНО: работаем ТОЛЬКО с копиями; рабочие prof и cfg остаются нетронутыми.
prof_ideal = zero_all_imu_errors(prof);
cfg_ideal  = zero_all_init_errors(cfg);

% Диагностику пишем на КАЖДОМ такте, чтобы корректно поймать максимумы
cfg_ideal.diag_decim = 1;

% Идеальные (нулевые) ошибки датчика: розыгрыш даст нули, но вызовы randn
% сохраняются, поэтому структура e_imu остаётся согласованной.
rng_ideal = RandStream('mt19937ar','Seed', cfg.seed);
e_ideal   = imu_draw_errors(prof_ideal, rng_ideal);
imu_ideal_meas = add_imu_errors(imu, prof_ideal, e_ideal, rng_ideal);

fprintf('\n=== ПРОВЕРКА ИДЕАЛЬНОСТИ ИЗМЕРЕНИЙ ===\n');
d_f = max(vecnorm(imu_ideal_meas.fb    - imu.fb,    2, 2));
d_w = max(vecnorm(imu_ideal_meas.wib_b - imu.wib_b, 2, 2));
fprintf('  макс |f_изм - f_ист|      = %.3e м/с²  (ожидаем ~0)\n', d_f);
fprintf('  макс |w_изм - w_ист|      = %.3e рад/с (ожидаем ~0)\n', d_w);
fprintf('  turn-on bias гиро/акс     = %.3e / %.3e\n', ...
        norm(e_ideal.gyro_turnon_bias), norm(e_ideal.accel_turnon_bias));

% =====================================================================
% Шаг T.4: TEST A — MECHANIZATION FLOOR (GNSS ОТКЛЮЧЁН ПОЛНОСТЬЮ)
% =====================================================================
% GNSS выключается через окно аутэйджа на ВЕСЬ полёт. Это не требует правки
% eskf_run: функция gnss_available вернёт false на каждом такте.
cfg_A = cfg_ideal;
cfg_A.outage = [ truth.t(1) - 1.0,  t_end + 1.0 ];

fprintf('\n=========================================================\n');
fprintf('  TEST A — MECHANIZATION FLOOR (без GNSS)\n');
fprintf('=========================================================\n');
tA = tic;
resA = eskf_run(imu_ideal_meas, truth, c, prof_ideal, cfg_A);
runtimeA = toc(tA);
mA = collect_floor_metrics(resA, runtimeA);
print_floor_metrics('TEST A (mechanization floor)', mA);

% Контроль: обновлений GNSS быть не должно
if isfield(resA, 'n_gnss')
    fprintf('  Обновлений GNSS: %d  [%s]\n', resA.n_gnss, ...
            ternary_str(resA.n_gnss == 0, 'OK — контур разомкнут', 'ОШИБКА: GNSS сработал'));
end

% =====================================================================
% Шаг T.5: TEST B — FULL LOOP IDEAL (идеальный GNSS + ESKF)
% =====================================================================
% Шум GNSS делаем очень малым, но НЕ нулевым: при sigma = 0 матрица R
% вырождается и инновационная ковариация S может стать плохо обусловленной.
cfg_B = cfg_ideal;
cfg_B.sigma_pos = 1e-3*ones(3,1);    % 1 мм
cfg_B.sigma_vel = 1e-5*ones(3,1);    % 0.01 мм/с

fprintf('\n=========================================================\n');
fprintf('  TEST B — FULL LOOP IDEAL (идеальный GNSS + ESKF)\n');
fprintf('=========================================================\n');
tB = tic;
resB = eskf_run(imu_ideal_meas, truth, c, prof_ideal, cfg_B);
runtimeB = toc(tB);
mB = collect_floor_metrics(resB, runtimeB);
print_floor_metrics('TEST B (full loop ideal)', mB);
if isfield(resB, 'n_gnss')
    fprintf('  Обновлений GNSS: %d\n', resB.n_gnss);
end

% =====================================================================
% Шаг T.6: СВОДКА
% =====================================================================
fprintf('\n=========================================================\n');
fprintf('  СВОДКА ЧИСЛЕННОГО ПОЛА  (f_nav = %d Гц, f_truth = %.0f Гц)\n', ...
        cfg.fnav, 1/truth.dt);
fprintf('=========================================================\n');
fprintf('  %-26s %14s %14s\n', 'Метрика', 'TEST A', 'TEST B');
fprintf('  %s\n', repmat('-',1,56));
fprintf('  %-26s %14.4e %14.4e\n', 'max |dr|, м',        mA.max_dr,   mB.max_dr);
fprintf('  %-26s %14.4e %14.4e\n', 'terminal |dr|, м',   mA.term_dr,  mB.term_dr);
fprintf('  %-26s %14.4e %14.4e\n', 'max |dv|, м/с',      mA.max_dv,   mB.max_dv);
fprintf('  %-26s %14.4e %14.4e\n', 'terminal |dv|, м/с', mA.term_dv,  mB.term_dv);
fprintf('  %-26s %14.4e %14.4e\n', 'max |dpsi|, мрад',   mA.max_dpsi, mB.max_dpsi);
fprintf('  %-26s %14.4e %14.4e\n', 'term |dpsi|, мрад',  mA.term_dpsi,mB.term_dpsi);
fprintf('  %-26s %14.2f %14.2f\n', 'runtime, с',         mA.runtime,  mB.runtime);
fprintf('=========================================================\n');

floor_test.cfg_fnav  = cfg.fnav;
floor_test.f_truth   = 1/truth.dt;
floor_test.testA     = mA;
floor_test.testB     = mB;
floor_test.resA      = resA;
floor_test.resB      = resB;

save('falco_numerical_floor.mat', 'floor_test');
fprintf('\nсохранено: falco_numerical_floor.mat\n');


% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% =========================================================================

function prof_out = zero_all_imu_errors(prof_in)
%ZERO_ALL_IMU_ERRORS  Копия профиля со всеми обнулёнными источниками ошибок.
%   Модель ошибок НЕ переписывается: обнуляются только СКО в копии структуры.
    prof_out = prof_in;

    % --- Гироскоп ---
    prof_out.gyro.turnon_bias_sigma  = 0;
    prof_out.gyro.inrun_bias_sigma   = 0;
    prof_out.gyro.arw                = 0;
    prof_out.gyro.scale_factor_sigma = 0;
    prof_out.gyro.misalign_axes      = 0;
    prof_out.gyro.misalign_ortho     = 0;
    prof_out.gyro.vrc                = 0;
    % Обнуляем ОШИБКУ КАЛИБРОВКИ g-sensitivity, но не сам коэффициент:
    % компенсация остаётся частью механизации и её остаток входит в пол.
    prof_out.gyro.g_sens_cal_sigma   = 0;

    % --- Акселерометр ---
    prof_out.accel.turnon_bias_sigma  = 0;
    prof_out.accel.inrun_bias_sigma   = 0;
    prof_out.accel.vrw                = 0;
    prof_out.accel.scale_factor_sigma = 0;
    prof_out.accel.misalign_axes      = 0;
    prof_out.accel.misalign_ortho     = 0;
    prof_out.accel.vrc                = 0;

    % --- Вибрация (страховка: при vrc=0 и так нет эффекта) ---
    prof_out.vib_grms = 0;
end


function cfg_out = zero_all_init_errors(cfg_in)
%ZERO_ALL_INIT_ERRORS  Копия конфигурации с нулевыми начальными ошибками
%   и нулевой аномалией гравитации.
    cfg_out = cfg_in;

    % Начальная ошибка выставки (детерминированная и статистическая формы)
    cfg_out.align_err   = zeros(3,1);
    cfg_out.align_sigma = zeros(3,1);

    % Аномалия гравитационного поля
    cfg_out.defl_vert_sigma = 0;
    cfg_out.grav_anom_sigma = 0;
    cfg_out.dg_model        = zeros(3,1);

    % ПРИМЕЧАНИЕ: P0 НЕ обнуляется. Это настройка фильтра, а не источник
    % ошибки. В TEST A обновлений нет вовсе, поэтому P на номинальную
    % траекторию не влияет. В TEST B заниженная P0 сделала бы фильтр
    % переуверенным, что исказило бы смысл теста.
end


function m = collect_floor_metrics(res, runtime)
%COLLECT_FLOOR_METRICS  Сбор метрик численного пола из результата прогона.
    dr_n   = vecnorm(res.dr,   2, 2);
    dv_n   = vecnorm(res.dv,   2, 2);
    dpsi_n = vecnorm(res.dpsi, 2, 2) * 1e3;   % в мрад

    m.max_dr    = max(dr_n);
    m.term_dr   = dr_n(end);
    m.max_dv    = max(dv_n);
    m.term_dv   = dv_n(end);
    m.max_dpsi  = max(dpsi_n);
    m.term_dpsi = dpsi_n(end);
    m.runtime   = runtime;

    [~, m.i_max_dr] = max(dr_n);
    m.t_max_dr      = res.t(m.i_max_dr);
    m.n_epochs      = numel(res.t);
end


function print_floor_metrics(title_str, m)
%PRINT_FLOOR_METRICS  Печать метрик одного теста.
    fprintf('  %s\n', title_str);
    fprintf('  %s\n', repmat('-',1,52));
    fprintf('    max |dr|          = %.6e м  (в t = %.2f с)\n', m.max_dr, m.t_max_dr);
    fprintf('    terminal |dr|     = %.6e м\n',   m.term_dr);
    fprintf('    max |dv|          = %.6e м/с\n', m.max_dv);
    fprintf('    terminal |dv|     = %.6e м/с\n', m.term_dv);
    fprintf('    max |dpsi|        = %.6e мрад\n', m.max_dpsi);
    fprintf('    terminal |dpsi|   = %.6e мрад\n', m.term_dpsi);
    fprintf('    эпох диагностики  = %d\n', m.n_epochs);
    fprintf('    runtime           = %.2f с\n', m.runtime);
end


function s = ternary_str(cond, a, b)
%TERNARY_STR  Выбор строки по условию.
    if cond, s = a; else, s = b; end
end


function print_time_grid(truth, fnav, f_gnss)
%PRINT_TIME_GRID  Проверка согласованности временных сеток (тест A1).
%
%   Печатает запрошенную и ФАКТИЧЕСКУЮ частоту навигатора. Расхождение
%   означает, что частота истины не кратна частоте навигатора и происходит
%   молчаливое округление числа отсчётов на такт.

    f_truth   = 1/truth.dt;
    ratio_raw = f_truth / fnav;
    step      = round(ratio_raw);
    dt_nav    = step * truth.dt;
    fnav_real = 1 / dt_nav;
    ratio_gn  = fnav / f_gnss;

    fprintf('\n=== ПРОВЕРКА ВРЕМЕННОЙ СЕТКИ (A1) ===\n');
    fprintf('  запрошенная f_nav   : %10.4f Гц\n', fnav);
    fprintf('  f_truth             : %10.4f Гц\n', f_truth);
    fprintf('  ratio f_truth/f_nav : %10.6f\n',    ratio_raw);
    fprintf('  step (отсчётов/такт): %10d\n',      step);
    fprintf('  dt_nav              : %10.6f с\n',  dt_nav);
    fprintf('  ФАКТИЧЕСКАЯ f_nav   : %10.4f Гц\n', fnav_real);
    fprintf('  ratio f_nav/f_gnss  : %10.4f\n',    ratio_gn);

    if abs(ratio_raw - step) > 1e-12
        fprintf('  СТАТУС: НЕКРАТНЫЕ ЧАСТОТЫ — фактическая f_nav НЕ равна запрошенной!\n');
        fprintf('          Запрошено %.3f Гц, реально %.3f Гц.\n', fnav, fnav_real);
    else
        fprintf('  СТАТУС: OK — частоты кратны, f_nav точная\n');
    end
    if abs(ratio_gn - round(ratio_gn)) > 1e-12
        fprintf('  ВНИМАНИЕ: f_nav/f_gnss не целое — такты GNSS не попадают\n');
        fprintf('            точно в такты навигатора.\n');
    end
end


function run_ned_checks(truth)
%RUN_NED_CHECKS  Регрессионные проверки NED-интерфейса.
%
%   Конвенция проекта:  v_e = C_e_ned * v_ned, то есть C_e_ned — матрица
%   NED -> ECEF, её столбцы суть орты [Север, Восток, Вниз] в ECEF.
%   Для величин, хранящихся ПОСТРОЧНО: ned_rows = ecef_rows * C_e_ned.

    fprintf('\n=== РЕГРЕССИОННЫЕ ПРОВЕРКИ NED ===\n');

    if ~isfield(truth,'C_e_ned0')
        fprintf('  truth.C_e_ned0 отсутствует (старый .mat) — проверки пропущены.\n');
        fprintf('  Перегенерируй falco_step7.mat через main_step7.\n');
        return;
    end

    C  = truth.C_e_ned0;
    Ct = truth.C_ned_e0;

    e_orth = norm(C'*C - eye(3));
    d_det  = det(C);
    e_tr   = norm(Ct - C');

    fprintf('  ||C_e_ned0''*C_e_ned0 - I|| = %.3e   [%s]\n', e_orth, ...
            pass_str(e_orth < 1e-12));
    fprintf('  det(C_e_ned0)              = %+.12f  [%s]\n', d_det, ...
            pass_str(abs(d_det - 1) < 1e-12));
    fprintf('  ||C_ned_e0 - C_e_ned0''||   = %.3e   [%s]\n', e_tr, ...
            pass_str(e_tr < 1e-15));

    if isfield(truth,'ned')
        e_start = norm(truth.ned(1,:));
        fprintf('  truth.ned(1,:)             = [%.3e %.3e %.3e]  [%s]\n', ...
                truth.ned(1,:), pass_str(e_start < 1e-6));
        if isfield(truth,'enu')
            % Связь: ned = [N E -U]
            e1 = max(abs(truth.ned(:,1) - truth.enu(:,2)));
            e2 = max(abs(truth.ned(:,2) - truth.enu(:,1)));
            e3 = max(abs(truth.ned(:,3) + truth.enu(:,3)));
            fprintf('  max|ned(:,1) - enu(:,2)|   = %.3e   [%s]\n', e1, pass_str(e1 < 1e-9));
            fprintf('  max|ned(:,2) - enu(:,1)|   = %.3e   [%s]\n', e2, pass_str(e2 < 1e-9));
            fprintf('  max|ned(:,3) + enu(:,3)|   = %.3e   [%s]\n', e3, pass_str(e3 < 1e-9));
        end
    end
end


function s = pass_str(ok)
%PASS_STR  Метка результата проверки.
    if ok, s = 'OK'; else, s = 'FAIL'; end
end
