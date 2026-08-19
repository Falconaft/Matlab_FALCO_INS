function plot_trajectory_v2(truth, imu, p)
%PLOT_TRAJECTORY_V2  Графики реалистичной траектории (Шаг 7).
%
%   Три фигуры:
%     1. Профиль полёта: высота, скорость, угол траектории, число Маха
%     2. Удельная сила и угловая скорость (что видит IMU) — КЛЮЧЕВОЕ
%     3. Траектория в вертикальной плоскости и в плане
%
%   Фазы подсвечены цветом: разгон / перелом тангажа / спуск / терминал.

    g0 = 9.80665;
    t  = truth.t;
    f_mag = vecnorm(imu.fb, 2, 2)/g0;
    w_deg = rad2deg(vecnorm(imu.wib_b, 2, 2));

    % =====================================================================
    % ФИГУРА 1: профиль полёта
    % =====================================================================
    figure('Color','w','Position',[40 40 1150 700]);

    ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, truth.alt/1000, 'LineWidth',1.7, 'Color',[0.12 0.31 0.61]);
    ylabel(ax,'Высота, км'); title(ax,'Профиль высоты','FontWeight','bold');
    shade_phases(ax, truth);

    ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, truth.Vmag, 'LineWidth',1.7, 'Color',[0.75 0.22 0.17]);
    ylabel(ax,'Скорость, м/с'); title(ax,'Модуль скорости','FontWeight','bold');
    shade_phases(ax, truth);

    ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, truth.gamma, 'LineWidth',1.7, 'Color',[0.49 0.24 0.60]);
    yline(ax, 0, ':', 'Color',[0.3 0.3 0.3]);
    ylabel(ax,'\gamma, град'); xlabel(ax,'Время, с');
    title(ax,'Угол наклона траектории','FontWeight','bold');
    shade_phases(ax, truth);

    ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, truth.Mach, 'LineWidth',1.7, 'Color',[0.15 0.55 0.30]);
    yline(ax, 1, '--', 'Color',[0.5 0.5 0.5], 'Label','M=1');
    ylabel(ax,'Число Маха'); xlabel(ax,'Время, с');
    title(ax,'Число Маха','FontWeight','bold');
    shade_phases(ax, truth);

    % =====================================================================
    % ФИГУРА 2: что видит IMU — КЛЮЧЕВОЕ для наблюдаемости
    % =====================================================================
    figure('Color','w','Position',[80 80 1150 700]);

    ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, f_mag, 'LineWidth',1.6, 'Color',[0.15 0.55 0.30]);
    set(ax,'YScale','log');
    yline(ax, 1, '--', 'Color',[0.5 0.5 0.5]);
    ylabel(ax,'|f|, g (лог. шкала)');
    title(ax,'Модуль удельной силы — определяет НАБЛЮДАЕМОСТЬ','FontWeight','bold');
    shade_phases(ax, truth);

    ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, imu.fb(:,1)/g0, 'LineWidth',1.4);
    plot(ax, t, imu.fb(:,2)/g0, 'LineWidth',1.4);
    plot(ax, t, imu.fb(:,3)/g0, 'LineWidth',1.4);
    ylabel(ax,'f^b, g'); legend(ax,{'X (вперёд)','Y (вправо)','Z (вниз)'},'Location','best');
    title(ax,'Удельная сила в осях тела','FontWeight','bold');
    shade_phases(ax, truth);

    ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, w_deg, 'LineWidth',1.6, 'Color',[0.85 0.45 0.10]);
    ylabel(ax,'|\omega|, °/с'); xlabel(ax,'Время, с');
    title(ax,'Модуль угловой скорости','FontWeight','bold');
    shade_phases(ax, truth);

    ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, truth.alpha, 'LineWidth',1.7, 'Color',[0.85 0.45 0.10]);
    ylabel(ax,'\alpha, град'); xlabel(ax,'Время, с');
    title(ax,'Программа угла атаки (сглаженная)','FontWeight','bold');
    shade_phases(ax, truth);

    % =====================================================================
    % ФИГУРА 3: геометрия траектории
    % =====================================================================
    figure('Color','w','Position',[120 120 1100 460]);

    ax = subplot(1,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    rng_km = hypot(truth.enu(:,1), truth.enu(:,2))/1000;
    plot(ax, rng_km, truth.alt/1000, 'LineWidth',1.9, 'Color',[0.12 0.31 0.61]);
    plot(ax, 0, truth.alt(1)/1000, 'go', 'MarkerFaceColor','g','MarkerSize',7);
    plot(ax, rng_km(end), truth.alt(end)/1000, 'rx', 'MarkerSize',11,'LineWidth',2);
    xlabel(ax,'Дальность, км'); ylabel(ax,'Высота, км');
    title(ax,'Траектория в вертикальной плоскости','FontWeight','bold');
    legend(ax,{'траектория','старт','удар'},'Location','best');

    ax = subplot(1,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    plot(ax, truth.enu(:,1)/1000, truth.enu(:,2)/1000, 'LineWidth',1.9, 'Color',[0.12 0.31 0.61]);
    plot(ax, 0, 0, 'go', 'MarkerFaceColor','g','MarkerSize',7);
    plot(ax, truth.enu(end,1)/1000, truth.enu(end,2)/1000, 'rx','MarkerSize',11,'LineWidth',2);
    xlabel(ax,'Восток, км'); ylabel(ax,'Север, км');
    title(ax,'Траектория в плане','FontWeight','bold');
end


function shade_phases(ax, truth)
%SHADE_PHASES  Подсветка фаз полёта БЕЗ искажения пределов осей.
%   Вызывается ПОСЛЕ построения данных: пределы фиксируются, патчи
%   рисуются внутри них, затем пределы восстанавливаются.

    cols = [0.82 0.91 1.00;    % фаза 1 — разгон (голубой)
            1.00 0.91 0.75;    % фаза 2 — перелом тангажа (оранжевый)
            0.85 0.94 0.85;    % фаза 3 — спуск (зелёный)
            1.00 0.85 0.85];   % фаза 4 — терминал (розовый)

    xl = xlim(ax);   yl = ylim(ax);

    for ph = 1:4
        m = truth.phase == ph;
        if ~any(m), continue; end
        i1 = find(m, 1);   i2 = find(m, 1, 'last');
        x1 = truth.t(i1);  x2 = truth.t(i2);
        pt = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], cols(ph,:), ...
                   'EdgeColor','none', 'HandleVisibility','off');
        uistack(pt, 'bottom');
    end

    xlim(ax, xl);   ylim(ax, yl);
    set(ax, 'Layer','top');
end
