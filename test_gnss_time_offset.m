%% TEST_GNSS_TIME_OFFSET  Проверка рассинхронизации GNSS/ИНС (B4)
%
%  МОДЕЛИРУЕТСЯ: некомпенсированная задержка измерения
%      t_meas = t_nav - cfg.gnss_time_offset
%  Измерение синтезируется по истине в момент t_meas, коррекция фильтра
%  выполняется в t_nav. Предсказание, H и update остаются в t_nav.
%
%  ЧАСТЬ 1 — ГЕОМЕТРИЧЕСКАЯ ПРОВЕРКА (без фильтра).
%      Смещение истинной антенны за время задержки должно составлять
%          |dr| ~ |v| * dt
%      Проверяется напрямую по truth, независимо от работы ESKF.
%
%  ЧАСТЬ 2 — SWEEP по сдвигам 0, 0.5, 1, 2, 5, 10 мс.
%      При одинаковом seed сравниваются terminal horizontal / 3D / velocity
%      ошибки и статистика невязок GNSS.
%      Ожидание: невязка позиции получает систематическую компоненту,
%      направленную ПРОТИВ вектора скорости, величиной ~|v|*dt. Это самый
%      чувствительный индикатор — виден раньше, чем эффект в конечной ошибке.
%
%  УСЛОВИЯ: B2 отключена, шум GNSS минимальный, идеальный IMU,
%  нулевые начальные ошибки, один и тот же seed во всех прогонах.
%
%  Требует falco_step7.mat.

clear; clc; close all;
load('falco_step7.mat', 'truth', 'imu', 'c', 'p');

% =====================================================================
% Шаг O.1: БАЗОВАЯ КОНФИГУРАЦИЯ
% =====================================================================
cfg = falco_config();
cfg.seed       = 5150;
cfg.f_diag     = 10;
cfg.diag_decim = round(cfg.fnav / cfg.f_diag);

% Шум GNSS минимальный, но НЕ нулевой: при sigma = 0 матрица R вырождается.
cfg.sigma_pos = 1e-3*ones(3,1);     % 1 мм
cfg.sigma_vel = 1e-5*ones(3,1);     % 0.01 мм/с

% Нулевые начальные ошибки
cfg.align_err   = zeros(3,1);
cfg.align_sigma = zeros(3,1);
cfg.P0_att      = [1e-3; 1e-3; 1e-3];

% B2 ОТКЛЮЧЕНА: изучаем эффект задержки в чистом виде
cfg.defl_vert_sigma = 0;
cfg.grav_anom_sigma = 0;
cfg.dg_model        = zeros(3,1);

t_end         = truth.t(end);
t_coast_start = t_end - p.coast_duration;
cfg.outage = [ 0.0            p.outage_boost_end;
               t_coast_start  t_end + 1.0 ];

% Идеальный датчик
prof_ideal = zero_imu_errors(imu_profile_pulse40_updated());

fprintf('=== TEST GNSS TIME OFFSET (B4) ===\n');
fprintf('Траектория: %.1f с, сетка истины %.0f Гц (шаг %.4f мс)\n', ...
        t_end, 1/truth.dt, truth.dt*1e3);
fprintf('Навигатор %.0f Гц, GNSS %.0f Гц, seed %d\n', ...
        cfg.fnav, cfg.f_gnss, cfg.seed);
fprintf('B2 отключена, шум GNSS %.0e м / %.0e м/с, IMU идеальный\n\n', ...
        cfg.sigma_pos(1), cfg.sigma_vel(1));

% =====================================================================
% Шаг O.2: ЧАСТЬ 1 — ГЕОМЕТРИЧЕСКАЯ ПРОВЕРКА ПОРЯДКА (без фильтра)
% =====================================================================
fprintf('%s\n', repmat('=',1,78));
fprintf('  ЧАСТЬ 1: геометрия — смещение ИСТИННОЙ АНТЕННЫ за время задержки\n');
fprintf('  Точка сравнения — фазовый центр антенны, а не центр IMU:\n');
fprintf('      r_ant(t) = truth.R(t) + Ceb(t)*lever\n');
fprintf('  Ожидаем |r_ant(t_nav) - r_ant(t_meas)| ~ |v|*dt\n');
fprintf('  Вклад плеча через поворот составляет |w x l|*dt и на порядки\n');
fprintf('  меньше основного члена, но учитывается корректно.\n');
fprintf('%s\n', repmat('=',1,78));

offsets_ms = [0, 0.5, 1, 2, 5, 10];
t_probe    = [20, 60, 120, 190];        % моменты с разной скоростью

% Позиция антенны в произвольный индекс сетки истины
ant_pos = @(idx) truth.R(idx,:)' + imu.Ceb(:,:,idx) * cfg.lever;

fprintf('  %8s', 'dt, мс');
for tp = t_probe
    fprintf(' | t=%3.0fс (V=%6.1f)', tp, interp1(truth.t, truth.Vmag, tp));
end
fprintf('\n  %s\n', repmat('-',1,8 + 22*numel(t_probe)));

for ms = offsets_ms
    n_sh = round(ms*1e-3 / truth.dt);
    fprintf('  %8.1f', ms);
    for tp = t_probe
        [~, ib] = min(abs(truth.t - tp));
        im = max(ib - n_sh, 1);
        % Сравниваем положения АНТЕННЫ, а не центра IMU
        d_act = norm(ant_pos(ib) - ant_pos(im));
        d_th  = truth.Vmag(ib) * ms*1e-3;
        fprintf(' | %8.4f/%7.4f м', d_act, d_th);
    end
    fprintf('\n');
end
fprintf('  (формат: фактическое смещение антенны / теоретическое |v|*dt)\n');

% Отдельно: насколько вклад плеча отличает антенну от центра IMU
fprintf('\n  Вклад ПОВОРОТА ПЛЕЧА (разность смещений антенны и центра IMU):\n');
fprintf('  %8s', 'dt, мс');
for tp = t_probe
    fprintf(' | t=%3.0fс', tp);
end
fprintf('\n  %s\n', repmat('-',1,8 + 11*numel(t_probe)));
for ms = offsets_ms(2:end)
    n_sh = round(ms*1e-3 / truth.dt);
    fprintf('  %8.1f', ms);
    for tp = t_probe
        [~, ib] = min(abs(truth.t - tp));
        im = max(ib - n_sh, 1);
        d_ant = norm(ant_pos(ib) - ant_pos(im));
        d_imu = norm(truth.R(ib,:) - truth.R(im,:));
        fprintf(' | %8.2e', abs(d_ant - d_imu));
    end
    fprintf('\n');
end
fprintf('\n');

% =====================================================================
% Шаг O.3: ЧАСТЬ 2 — SWEEP ПО СДВИГАМ
% =====================================================================
fprintf('%s\n', repmat('=',1,72));
fprintf('  ЧАСТЬ 2: sweep по временным сдвигам (одинаковый seed)\n');
fprintf('%s\n', repmat('=',1,72));

n_off = numel(offsets_ms);
sw = struct('offset_ms',cell(n_off,1), 'err_horiz',[], 'err_3d',[], ...
            'err_v',[], 'innov_pos_mean',[], 'innov_pos_std',[], ...
            'innov_vel_mean',[], 'innov_vel_std',[], 'mismatch_med',[], ...
            'n_gnss',[], 'ba_est_norm',[], ...
            'pre_horiz',[], 'pre_3d',[], 'pre_dv',[], 'pre_dpsi',[], ...
            'pre_ba',[], 'zpar_mean',[], 'zpar_frac_neg',[], ...
            'zpar_mean_early',[], 'zpar_mean_steady',[]);

for i = 1:n_off
    ms = offsets_ms(i);

    cfg_i = cfg;
    cfg_i.gnss_time_offset_enable = (ms > 0);
    cfg_i.gnss_time_offset        = ms * 1e-3;

    % Один и тот же seed -> одинаковая последовательность шума
    rng_i    = RandStream('mt19937ar','Seed', cfg.seed);
    e_i      = imu_draw_errors(prof_ideal, rng_i);
    imu_meas = add_imu_errors(imu, prof_ideal, e_i, rng_i);

    res = eskf_run(imu_meas, truth, c, prof_ideal, cfg_i);

    sw(i).offset_ms      = ms;
    sw(i).err_horiz      = res.err_horiz_final;
    sw(i).err_3d         = res.err_3d_final;
    sw(i).err_v          = norm(res.dv(end,:));
    sw(i).n_gnss         = res.n_gnss;
    sw(i).innov_pos_mean = mean(vecnorm(res.gnss_innov(:,1:3),2,2));
    sw(i).innov_pos_std  = std(vecnorm(res.gnss_innov(:,1:3),2,2));
    sw(i).innov_vel_mean = mean(vecnorm(res.gnss_innov(:,4:6),2,2));
    sw(i).innov_vel_std  = std(vecnorm(res.gnss_innov(:,4:6),2,2));
    if isfield(res,'gnss_mismatch') && ~isempty(res.gnss_mismatch)
        sw(i).mismatch_med = median(vecnorm(res.gnss_mismatch,2,2));
    else
        sw(i).mismatch_med = NaN;
    end

    % ФИКТИВНЫЙ BIAS АКСЕЛЕРОМЕТРА.
    % Датчик идеальный, истинный bias РОВНО НОЛЬ. Поэтому любая ненулевая
    % оценка ba целиком порождена фильтром: он пытается объяснить
    % систематическую невязку от временного сдвига через bias. Рост этой
    % величины по dt означает, что часть эффекта задержки поглощается
    % в неверную оценку bias, а не проявляется в конечной ошибке.
    sw(i).ba_est_norm = norm(res.nav.ba);

    % ---- СОСТОЯНИЕ НЕПОСРЕДСТВЕННО ПЕРЕД ТЕРМИНАЛЬНЫМ КОАСТОМ ----
    % Ключевая точка цепочки: time offset -> ошибка перед коастом ->
    % ошибка после 60 с чистой ИНС. Позволяет разделить, сколько ошибки
    % внесено самой задержкой, а сколько накоплено на бесспутниковом участке.
    [~, k_pre]      = min(abs(res.t - t_coast_start));
    sw(i).pre_horiz = hypot(res.dr_ned(k_pre,1), res.dr_ned(k_pre,2));
    sw(i).pre_3d    = norm(res.dr(k_pre,:));
    sw(i).pre_dv    = norm(res.dv(k_pre,:));
    sw(i).pre_dpsi  = norm(res.dpsi(k_pre,:)) * 1e3;      % мрад
    sw(i).pre_ba    = norm(res.ba_est(k_pre,:));

    % ---- ЗНАКОВАЯ ПРОЕКЦИЯ НЕВЯЗКИ НА НАПРАВЛЕНИЕ ИСТИННОЙ СКОРОСТИ ----
    % z_par = <z_r, v_hat>. Наивное ожидание для положительного сдвига:
    % z_par < 0, поскольку измерение относится к прошлому и лежит "позади".
    %
    % ОДНАКО в установившемся режиме первый порядок СОКРАЩАЕТСЯ: фильтр
    % уже отслеживает запаздывающую истину, механизация ведёт состояние
    % вперёд истинным ускорением, и члены -v*dt и +v*dt взаимно гасятся,
    % оставляя лишь второй порядок (изменение ускорения за такт, усиление
    % Калмана < 1). Поэтому знак информативен, но не предопределён.
    % Наиболее показателен участок СРАЗУ ПОСЛЕ возобновления GNSS, где
    % установившийся режим ещё не наступил — он выводится отдельно.
    v_at_gnss = interp1(truth.t, truth.V, res.gnss_t, 'linear');
    v_norm    = vecnorm(v_at_gnss, 2, 2);
    v_hat     = v_at_gnss ./ max(v_norm, eps);
    z_r       = res.gnss_innov(:,1:3);
    z_par     = sum(z_r .* v_hat, 2);

    sw(i).zpar_mean     = mean(z_par);
    sw(i).zpar_frac_neg = mean(z_par < 0);

    % Разделение: первые 5 с после возобновления GNSS (переходный режим)
    % против остального окна (установившийся)
    t_gnss_on = p.outage_boost_end;
    i_early   = res.gnss_t <= t_gnss_on + 5.0;
    i_steady  = ~i_early;
    if any(i_early),  sw(i).zpar_mean_early  = mean(z_par(i_early));
    else,             sw(i).zpar_mean_early  = NaN; end
    if any(i_steady), sw(i).zpar_mean_steady = mean(z_par(i_steady));
    else,             sw(i).zpar_mean_steady = NaN; end

    fprintf('  сдвиг %5.1f мс: n_shift=%2d, обновлений %4d, гориз. %8.4f м\n', ...
            ms, res.gnss_n_shift, res.n_gnss, res.err_horiz_final);

    if ms == 0
        res_ref = res;
    end
end

% =====================================================================
% Шаг O.4: СВОДНАЯ ТАБЛИЦА
% =====================================================================
fprintf('\n%s\n', repmat('=',1,100));
fprintf('  СВОДКА SWEEP\n');
fprintf('%s\n', repmat('=',1,100));
fprintf('  %7s %11s %11s %12s %13s %13s %11s %13s\n', ...
        'dt, мс', 'гориз., м', '3D, м', '|dv|, м/с', 'невязка r, м', ...
        'невязка v, м/с', 'рассогл., м', '|ba| оц., мг');
fprintf('  %s\n', repmat('-',1,100));
g0_loc = 9.80665;
for i = 1:n_off
    fprintf('  %7.1f %11.4f %11.4f %12.5f %13.5f %13.6f %11.4f %13.5f\n', ...
            sw(i).offset_ms, sw(i).err_horiz, sw(i).err_3d, sw(i).err_v, ...
            sw(i).innov_pos_mean, sw(i).innov_vel_mean, sw(i).mismatch_med, ...
            sw(i).ba_est_norm/g0_loc*1e3);
end
fprintf('%s\n', repmat('=',1,100));
fprintf('  невязка r / v — СРЕДНИЙ модуль невязки GNSS по всем обновлениям\n');
fprintf('  рассогл. — медиана |r_ant_true(t_meas) - r_ant_true(t_nav)|\n');
fprintf('  |ba| оц. — оценка bias акселерометра; ИСТИННЫЙ bias здесь РОВНО 0,\n');
fprintf('             поэтому любое ненулевое значение = поглощение сдвига\n');
fprintf('             фильтром в фиктивный bias\n');

% --- Состояние перед терминальным коастом ---
fprintf('\n%s\n', repmat('=',1,96));
fprintf('  СОСТОЯНИЕ ПЕРЕД ТЕРМИНАЛЬНЫМ КОАСТОМ (t = %.1f с) и ПОСЛЕ НЕГО\n', t_coast_start);
fprintf('%s\n', repmat('=',1,96));
fprintf('  %7s | %10s %10s %11s %11s %11s | %10s %10s\n', ...
        'dt, мс', 'гориз.,м', '3D, м', '|dv|,м/с', '|dpsi|,мрад', '|ba|,мг', ...
        'гориз.кон', 'рост, м');
fprintf('  %s\n', repmat('-',1,92));
for i = 1:n_off
    fprintf('  %7.1f | %10.4f %10.4f %11.6f %11.5f %11.6f | %10.4f %10.4f\n', ...
            sw(i).offset_ms, sw(i).pre_horiz, sw(i).pre_3d, sw(i).pre_dv, ...
            sw(i).pre_dpsi, sw(i).pre_ba/9.80665*1e3, ...
            sw(i).err_horiz, sw(i).err_horiz - sw(i).pre_horiz);
end
fprintf('%s\n', repmat('=',1,96));
fprintf('  рост = насколько ошибка выросла за 60 с чистой ИНС\n');

% --- Знаковая проекция невязки на направление скорости ---
fprintf('\n%s\n', repmat('=',1,80));
fprintf('  ЗНАКОВАЯ ПРОЕКЦИЯ НЕВЯЗКИ НА НАПРАВЛЕНИЕ ИСТИННОЙ СКОРОСТИ\n');
fprintf('  z_par = <z_r, v_hat>;  ожидание для dt > 0: измерение "позади"\n');
fprintf('%s\n', repmat('=',1,80));
fprintf('  %7s %14s %14s %16s %12s\n', ...
        'dt, мс', 'среднее, м', 'перех. (5с), м', 'установ., м', 'доля z<0');
fprintf('  %s\n', repmat('-',1,76));
for i = 1:n_off
    fprintf('  %7.1f %14.6f %14.6f %16.6f %11.1f%%\n', ...
            sw(i).offset_ms, sw(i).zpar_mean, sw(i).zpar_mean_early, ...
            sw(i).zpar_mean_steady, 100*sw(i).zpar_frac_neg);
end
fprintf('%s\n', repmat('=',1,80));
fprintf('  перех. — первые 5 с после возобновления GNSS (t > %.1f с)\n', ...
        p.outage_boost_end);
fprintf('  установ. — остальное окно GNSS\n');
fprintf('  ПРИМЕЧАНИЕ: в установившемся режиме первый порядок сокращается\n');
fprintf('  (фильтр уже отслеживает запаздывающую истину), поэтому знак там\n');
fprintf('  определяется вторым порядком и может не совпасть с наивным.\n');

% --- Рост относительно нулевого сдвига ---
fprintf('\n  РОСТ ОТНОСИТЕЛЬНО dt = 0:\n');
fprintf('  %7s %14s %14s %16s\n', 'dt, мс', 'гориз., м', 'прирост, м', 'невязка r, м');
fprintf('  %s\n', repmat('-',1,56));
for i = 1:n_off
    fprintf('  %7.1f %14.4f %14.4f %16.5f\n', sw(i).offset_ms, sw(i).err_horiz, ...
            sw(i).err_horiz - sw(1).err_horiz, sw(i).innov_pos_mean);
end

% =====================================================================
% Шаг O.5: ГРАФИКИ
% =====================================================================
off = [sw.offset_ms]';
figure('Color','w','Position',[60 60 1050 640]);

ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, off, [sw.err_horiz]', '-o', 'LineWidth',1.7, 'MarkerFaceColor','auto');
plot(ax, off, [sw.err_3d]',    '-s', 'LineWidth',1.4);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'Ошибка на ударе, м');
title(ax,'Конечная ошибка от временного сдвига','FontWeight','bold');
legend(ax,{'горизонтальная','3D'},'Location','best');

ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, off, [sw.innov_pos_mean]', '-o', 'LineWidth',1.7);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'Средняя |невязка позиции|, м');
title(ax,'Невязка GNSS — самый чувствительный индикатор','FontWeight','bold');

ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, off, [sw.innov_vel_mean]', '-o', 'LineWidth',1.7, 'Color',[0.8 0.3 0.2]);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'Средняя |невязка скорости|, м/с');
title(ax,'Невязка по скорости','FontWeight','bold');

ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, off, [sw.mismatch_med]', '-o', 'LineWidth',1.7, 'Color',[0.3 0.6 0.3]);
xlabel(ax,'Сдвиг GNSS, мс'); ylabel(ax,'Медиана рассогласования, м');
title(ax,'Геометрическое рассогласование истины','FontWeight','bold');

save('falco_test_gnss_offset.mat', 'sw', 'offsets_ms', 'cfg');
fprintf('\nсохранено: falco_test_gnss_offset.mat\n');


% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% =========================================================================

function prof_out = zero_imu_errors(prof_in)
%ZERO_IMU_ERRORS  Копия профиля со всеми обнулёнными источниками ошибок.
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
