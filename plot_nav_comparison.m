function plot_nav_comparison(truth, res, p, cfg, c)
%PLOT_NAV_COMPARISON  Сравнение истинной и навигационной траекторий.
%
%   Две фигуры:
%     Фигура 1 — геометрия:
%        (а) 3D truth vs nav в NED
%        (б) наземный трек Восток/Север
%        (в) профиль высоты по времени
%        (г) отклонение nav от truth в плане (увеличенно)
%     Фигура 2 — ошибки:
%        (а) ошибка позиции по компонентам NED
%        (б) горизонтальная ошибка по времени
%        (в) вертикальная ошибка по времени
%        (г) ошибка скорости по компонентам NED
%
%   Все графики строятся в NED (внешний интерфейс). Внутренний
%   вычислительный core остаётся в ECEF — здесь только преобразование вывода.
%
%   Отмечаются участки: работа двигателя, аутэйдж GNSS на разгоне,
%   терминальный коаст.
%
%   Вход:
%     truth - истинная траектория (нужны .t, .R, .r0, .ned или .Cen)
%     res   - результат eskf_run (нужны .t, .r_nav, .dr_ned, .dv_ned, .outage)
%     p     - параметры траектории (нужны .t_burn, .coast_duration)
%     cfg   - конфигурация (опционально, для подписи)

    if nargin < 4, cfg = struct(); end

    % =====================================================================
    % Шаг P.1: ПОДГОТОВКА ДАННЫХ В NED
    % =====================================================================
    % Матрица NED в точке старта: из truth либо из res (там же лежит копия),
    % либо выводится из Cen перестановкой столбцов (обратная совместимость).
    if isfield(truth, 'C_e_ned0')
        C_e_ned0 = truth.C_e_ned0;
    elseif isfield(res, 'C_e_ned0')
        C_e_ned0 = res.C_e_ned0;
    else
        C_e_ned0 = [truth.Cen(:,2), truth.Cen(:,1), -truth.Cen(:,3)];
    end

    r0 = truth.r0(:)';

    % Истинная траектория в NED (на СВОЕЙ временной сетке)
    if isfield(truth, 'ned')
        ned_truth = truth.ned;
    else
        ned_truth = (truth.R - r0) * C_e_ned0;
    end

    % Навигационная траектория в NED (на сетке диагностики)
    ned_nav = (res.r_nav - r0) * C_e_ned0;
    alt_nav = zeros(size(res.t));

    for k = 1:numel(res.t)
        [~, ~, alt_nav(k)] = ecef2lla(res.r_nav(k,:)', c);
    end
    
    
    % Истинная траектория, интерполированная на сетку диагностики —
    % нужна, чтобы рисовать обе кривые в одних точках
    ned_truth_i = interp1(truth.t, ned_truth, res.t, 'linear');
    alt_truth_i = interp1(truth.t, truth.alt, res.t, 'linear');
    
    t   = res.t;
    N_t = ned_truth_i(:,1);   E_t = ned_truth_i(:,2);   D_t = ned_truth_i(:,3);
    N_n = ned_nav(:,1);       E_n = ned_nav(:,2);       D_n = ned_nav(:,3);

    err_h = hypot(res.dr_ned(:,1), res.dr_ned(:,2));
    err_v = res.dr_ned(:,3);

    % Границы участков
    t_burn        = get_field(p, 't_burn', NaN);
    t_coast_start = t(end) - get_field(p, 'coast_duration', 0);

    % =====================================================================
    % ФИГУРА 1: ГЕОМЕТРИЯ
    % =====================================================================
    figure('Color','w','Position',[40 40 1150 720]);

    % --- (а) 3D ---
    ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot3(ax, E_t/1000, N_t/1000, -D_t/1000, 'LineWidth',2.0, 'Color',[0.15 0.35 0.70]);
    plot3(ax, E_n/1000, N_n/1000, -D_n/1000, '--', 'LineWidth',1.4, 'Color',[0.85 0.25 0.15]);
    plot3(ax, E_t(1)/1000, N_t(1)/1000, -D_t(1)/1000, 'go','MarkerFaceColor','g','MarkerSize',7);
    plot3(ax, E_t(end)/1000, N_t(end)/1000, -D_t(end)/1000, 'kx','MarkerSize',11,'LineWidth',2);
    xlabel(ax,'Восток, км'); ylabel(ax,'Север, км'); zlabel(ax,'Высота, км');
    title(ax,'Траектория 3D (NED): истина и навигация','FontWeight','bold');
    legend(ax,{'истина','навигация','старт','удар'},'Location','best');
    view(ax, -35, 22);

    % --- (б) наземный трек ---
    ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    plot(ax, E_t/1000, N_t/1000, 'LineWidth',2.0, 'Color',[0.15 0.35 0.70]);
    plot(ax, E_n/1000, N_n/1000, '--', 'LineWidth',1.4, 'Color',[0.85 0.25 0.15]);
    plot(ax, 0, 0, 'go','MarkerFaceColor','g','MarkerSize',7);
    plot(ax, E_t(end)/1000, N_t(end)/1000, 'kx','MarkerSize',11,'LineWidth',2);
    xlabel(ax,'Восток, км'); ylabel(ax,'Север, км');
    title(ax,'Наземный трек','FontWeight','bold');
    legend(ax,{'истина','навигация'},'Location','best');

    % --- (в) геодезическая высота ---
    ax = subplot(2,2,3);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    
    plot(ax, t, alt_truth_i/1000, 'LineWidth',2.0);
    plot(ax, t, alt_nav/1000, '--', 'LineWidth',1.4);
    
    xlabel(ax,'Время, с');
    ylabel(ax,'Высота, км');
    title(ax,'Геодезическая высота: истина и INS/GNSS','FontWeight','bold');
    
    mark_phases(ax, res, t_burn, t_coast_start);
    legend(ax,{'истина','INS/GNSS'},'Location','best');

    % --- (г) отклонение в плане ---
    ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    plot(ax, res.dr_ned(:,2), res.dr_ned(:,1), 'LineWidth',1.3, 'Color',[0.55 0.20 0.60]);
    plot(ax, 0, 0, 'k+','MarkerSize',13,'LineWidth',1.6);
    plot(ax, res.dr_ned(end,2), res.dr_ned(end,1), 'ro','MarkerFaceColor','r','MarkerSize',7);
    xlabel(ax,'Ошибка Восток, м'); ylabel(ax,'Ошибка Север, м');
    title(ax, sprintf('Отклонение навигации в плане (финал %.2f м)', err_h(end)), ...
          'FontWeight','bold');
    legend(ax,{'траектория ошибки','истина','финал'},'Location','best');

    % =====================================================================
    % ФИГУРА 2: ОШИБКИ
    % =====================================================================
    figure('Color','w','Position',[90 90 1150 720]);
    lbl = {'Север','Восток','Вниз'};
    col = lines(3);

    % --- (а) ошибка позиции по компонентам ---
    ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for i = 1:3
        plot(ax, t, res.dr_ned(:,i), 'LineWidth',1.3, 'Color',col(i,:));
    end
    ylabel(ax,'Ошибка позиции, м'); xlabel(ax,'Время, с');
    title(ax,'Ошибка позиции по осям NED','FontWeight','bold');
    mark_phases(ax, res, t_burn, t_coast_start);
    legend(ax, lbl, 'Location','best');

    % --- (б) горизонтальная ошибка ---
    ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, err_h, 'LineWidth',1.6, 'Color',[0.75 0.20 0.15]);
    ylabel(ax,'Горизонтальная ошибка, м'); xlabel(ax,'Время, с');
    title(ax, sprintf('Горизонтальная ошибка (финал %.3f м, макс %.3f м)', ...
          err_h(end), max(err_h)), 'FontWeight','bold');
    mark_phases(ax, res, t_burn, t_coast_start);

    % --- (в) вертикальная ошибка ---
    ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, err_v, 'LineWidth',1.6, 'Color',[0.20 0.50 0.30]);
    yline(ax, 0, ':', 'Color',[0.4 0.4 0.4]);
    ylabel(ax,'Ошибка по вертикали (Вниз), м'); xlabel(ax,'Время, с');
    title(ax, sprintf('Вертикальная ошибка (финал %.3f м)', err_v(end)), ...
          'FontWeight','bold');
    mark_phases(ax, res, t_burn, t_coast_start);

    % --- (г) ошибка скорости ---
    ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(res,'dv_ned')
        for i = 1:3
            plot(ax, t, res.dv_ned(:,i), 'LineWidth',1.3, 'Color',col(i,:));
        end
        legend(ax, lbl, 'Location','best');
    else
        plot(ax, t, vecnorm(res.dv,2,2), 'LineWidth',1.4);
        legend(ax, {'|dv| (ECEF)'}, 'Location','best');
    end
    ylabel(ax,'Ошибка скорости, м/с'); xlabel(ax,'Время, с');
    title(ax,'Ошибка скорости по осям NED','FontWeight','bold');
    mark_phases(ax, res, t_burn, t_coast_start);
end


% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% =========================================================================

function mark_phases(ax, res, t_burn, t_coast_start)
%MARK_PHASES  Подсветка участков БЕЗ искажения пределов осей.
%   Вызывается ПОСЛЕ построения данных: пределы фиксируются, патчи рисуются
%   внутри них, затем пределы восстанавливаются. Иначе патч растягивает
%   ось Y и данные схлопываются.

    xl = xlim(ax);   yl = ylim(ax);

    % --- Окна отсутствия GNSS (розовые) ---
    if isfield(res, 'outage')
        for i = 1:size(res.outage,1)
            x1 = max(res.outage(i,1), xl(1));
            x2 = min(res.outage(i,2), xl(2));
            if x2 > x1
                pt = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
                           [1 0.86 0.89], 'EdgeColor','none','HandleVisibility','off');
                uistack(pt,'bottom');
            end
        end
    end

    % --- Работа двигателя (голубая) ---
    if ~isnan(t_burn) && t_burn > xl(1)
        x2 = min(t_burn, xl(2));
        pt = patch(ax, [xl(1) x2 x2 xl(1)], [yl(1) yl(1) yl(2) yl(2)], ...
                   [0.85 0.92 1.0], 'EdgeColor','none','HandleVisibility','off');
        uistack(pt,'bottom');
    end

    % --- Границы участков (линии) ---
    if ~isnan(t_burn)
        xline(ax, t_burn, '-', 'Color',[0.2 0.4 0.8], 'LineWidth',1.1, ...
              'HandleVisibility','off');
    end
    if t_coast_start > xl(1) && t_coast_start < xl(2)
        xline(ax, t_coast_start, '-', 'Color',[0.8 0.2 0.2], 'LineWidth',1.1, ...
              'HandleVisibility','off');
    end

    xlim(ax, xl);   ylim(ax, yl);
    set(ax,'Layer','top');
end


function v = get_field(s, name, default)
%GET_FIELD  Чтение поля структуры со значением по умолчанию.
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end
