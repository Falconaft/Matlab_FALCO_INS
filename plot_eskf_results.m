function plot_eskf_results(res, truth, ba_true, bg_true)
%PLOT_ESKF_RESULTS  Диагностические графики работы ESKF.
%
%   Шесть фигур для визуального контроля сходимости фильтра:
%     1. Ошибка позиции (ENU) с коридором +-3*sigma и окнами аутэйджа
%     2. Ошибка скорости (ECEF) с коридором +-3*sigma
%     3. Ошибка ориентации с коридором +-3*sigma
%     4. Сходимость оценок bias (оценка vs истина)
%     5. Эволюция СКО всех состояний (log) - "обучение" фильтра
%     6. Невязки GNSS с границами +-3*sigma (проверка состоятельности)
%
%   ВАЖНО о подсветке аутэйджа: прямоугольники рисуются ПОСЛЕ данных, и
%   пределы осей восстанавливаются вручную. Иначе патч с большой высотой
%   растягивает ось Y и данные схлопываются в линию на нуле.
%
%   Вход:
%     res      - результат eskf_run
%     truth    - истинная траектория (для Cen)
%     ba_true  - истинный bias акселерометра (3x1), опц. []
%     bg_true  - истинный bias гироскопа (3x1), опц. []

    if nargin < 3, ba_true = []; end
    if nargin < 4, bg_true = []; end

    t    = res.t;
    col  = lines(3);
    g0   = 9.80665;
    deg  = pi/180;
    lbl_enu = {'Восток','Север','Верх'};
    lbl_xyz = {'X','Y','Z'};

    % Перенос СКО позиции из ECEF в ENU (по диагонали, приближённо)
    sig_enu = sqrt((res.sig_r.^2) * (truth.Cen.^2));

    % =====================================================================
    % ФИГУРА 1: ошибка позиции в ENU
    % =====================================================================
    figure('Color','w','Position',[50 50 1000 620]);
    for i = 1:3
        ax = subplot(3,1,i); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        plot(ax, t,  3*sig_enu(:,i), '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1);
        plot(ax, t, -3*sig_enu(:,i), '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1);
        plot(ax, t, res.dr_enu(:,i), 'LineWidth',1.5, 'Color',col(i,:));
        ylabel(ax, sprintf('%s, м', lbl_enu{i}));
        shade_outage(ax, res.outage);       % подсветка ПОСЛЕ данных
        if i==1
            title(ax, ['Ошибка позиции (ENU): линия — ошибка, пунктир — \pm3\sigma, ' ...
                       'розовым — нет GNSS'], 'FontWeight','bold');
        end
        if i==3, xlabel(ax,'Время, с'); end
    end

    % =====================================================================
    % ФИГУРА 2: ошибка скорости
    % =====================================================================
    figure('Color','w','Position',[80 80 1000 620]);
    for i = 1:3
        ax = subplot(3,1,i); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        plot(ax, t,  3*res.sig_v(:,i), '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1);
        plot(ax, t, -3*res.sig_v(:,i), '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1);
        plot(ax, t, res.dv(:,i), 'LineWidth',1.5, 'Color',col(i,:));
        ylabel(ax, sprintf('V%s, м/с', lbl_xyz{i}));
        shade_outage(ax, res.outage);
        if i==1, title(ax,'Ошибка скорости (ECEF) и коридор \pm3\sigma','FontWeight','bold'); end
        if i==3, xlabel(ax,'Время, с'); end
    end

    % =====================================================================
    % ФИГУРА 3: ошибка ориентации
    % =====================================================================
    figure('Color','w','Position',[110 110 1000 620]);
    for i = 1:3
        ax = subplot(3,1,i); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        plot(ax, t,  3*res.sig_p(:,i)*1e3, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1);
        plot(ax, t, -3*res.sig_p(:,i)*1e3, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1);
        plot(ax, t, res.dpsi(:,i)*1e3, 'LineWidth',1.5, 'Color',col(i,:));
        ylabel(ax, sprintf('\\delta\\psi_%s, мрад', lbl_xyz{i}));
        shade_outage(ax, res.outage);
        if i==1, title(ax,'Ошибка ориентации и коридор \pm3\sigma','FontWeight','bold'); end
        if i==3, xlabel(ax,'Время, с'); end
    end

    % =====================================================================
    % ФИГУРА 4: сходимость оценок bias
    % =====================================================================
    figure('Color','w','Position',[140 140 1050 700]);
    for i = 1:3
        % --- акселерометр ---
        ax = subplot(3,2,2*i-1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        plot(ax, t, (res.ba_est(:,i)+3*res.sig_ba(:,i))/g0*1e3, ':', 'Color',[0.6 0.6 0.6]);
        plot(ax, t, (res.ba_est(:,i)-3*res.sig_ba(:,i))/g0*1e3, ':', 'Color',[0.6 0.6 0.6]);
        plot(ax, t, res.ba_est(:,i)/g0*1e3, 'LineWidth',1.6, 'Color',col(1,:));
        if ~isempty(ba_true)
            yline(ax, ba_true(i)/g0*1e3, 'r--', 'LineWidth',1.4);
        end
        ylabel(ax, sprintf('b_{a%s}, мг', lbl_xyz{i}));
        shade_outage(ax, res.outage);
        if i==1, title(ax,'Bias акселерометра: оценка и \pm3\sigma (красным — истина)','FontWeight','bold'); end
        if i==3, xlabel(ax,'Время, с'); end

        % --- гироскоп ---
        ax = subplot(3,2,2*i); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        plot(ax, t, (res.bg_est(:,i)+3*res.sig_bg(:,i))/deg*3600, ':', 'Color',[0.6 0.6 0.6]);
        plot(ax, t, (res.bg_est(:,i)-3*res.sig_bg(:,i))/deg*3600, ':', 'Color',[0.6 0.6 0.6]);
        plot(ax, t, res.bg_est(:,i)/deg*3600, 'LineWidth',1.6, 'Color',col(2,:));
        if ~isempty(bg_true)
            yline(ax, bg_true(i)/deg*3600, 'r--', 'LineWidth',1.4);
        end
        ylabel(ax, sprintf('b_{g%s}, °/ч', lbl_xyz{i}));
        shade_outage(ax, res.outage);
        if i==1, title(ax,'Bias гироскопа: оценка и \pm3\sigma (красным — истина)','FontWeight','bold'); end
        if i==3, xlabel(ax,'Время, с'); end
    end

    % =====================================================================
    % ФИГУРА 5: эволюция СКО (обучение фильтра), лог. шкала
    % =====================================================================
    figure('Color','w','Position',[170 170 1000 640]);

    ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, res.sig_r, 'LineWidth',1.3); set(ax,'YScale','log');
    ylabel(ax,'\sigma позиции, м'); title(ax,'Неопределённость позиции','FontWeight','bold');
    legend(ax, lbl_xyz, 'Location','best');
    shade_outage(ax, res.outage);

    ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, res.sig_v, 'LineWidth',1.3); set(ax,'YScale','log');
    ylabel(ax,'\sigma скорости, м/с'); title(ax,'Неопределённость скорости','FontWeight','bold');
    shade_outage(ax, res.outage);

    ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, res.sig_p*1e3, 'LineWidth',1.3); set(ax,'YScale','log');
    ylabel(ax,'\sigma ориентации, мрад'); xlabel(ax,'Время, с');
    title(ax,'Неопределённость ориентации','FontWeight','bold');
    shade_outage(ax, res.outage);

    ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, res.sig_ba/g0*1e3, 'LineWidth',1.3); set(ax,'YScale','log');
    ylabel(ax,'\sigma bias акс., мг'); xlabel(ax,'Время, с');
    title(ax,'Неопределённость bias акселерометра','FontWeight','bold');
    shade_outage(ax, res.outage);

    % =====================================================================
    % ФИГУРА 6: невязки GNSS (проверка состоятельности фильтра)
    % =====================================================================
    figure('Color','w','Position',[200 200 1000 560]);
    if isfield(res,'gnss_t') && ~isempty(res.gnss_t)
        tg = res.gnss_t;
        ax = subplot(2,1,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        h = plot(ax, tg, res.gnss_innov(:,1:3), '.-', 'MarkerSize',10, 'LineWidth',0.6);
        ylabel(ax,'Невязка позиции, м'); legend(ax, h, lbl_xyz, 'Location','best');
        title(ax, sprintf(['Невязки GNSS (%d обновлений): должны лежать в пределах шума ' ...
                           'и не иметь тренда'], res.n_gnss), 'FontWeight','bold');
        shade_outage(ax, res.outage);

        ax = subplot(2,1,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        h = plot(ax, tg, res.gnss_innov(:,4:6), '.-', 'MarkerSize',10, 'LineWidth',0.6);
        ylabel(ax,'Невязка скорости, м/с'); xlabel(ax,'Время, с');
        legend(ax, h, lbl_xyz, 'Location','best');
        shade_outage(ax, res.outage);
    else
        text(0.5,0.5,'Нет записанных невязок GNSS','HorizontalAlignment','center');
        axis off;
    end
end


function shade_outage(ax, outage)
%SHADE_OUTAGE  Подсветка окон отсутствия GNSS БЕЗ искажения пределов осей.
%
%   Ключевой момент: вызывается ПОСЛЕ построения данных. Пределы осей
%   фиксируются, патч рисуется строго внутри них, затем пределы
%   восстанавливаются. Иначе патч растягивает ось Y (был баг: высота
%   патча 1e6 схлопывала реальные данные в линию на нуле).

    xl = xlim(ax);
    yl = ylim(ax);

    for i = 1:size(outage,1)
        x1 = max(outage(i,1), xl(1));
        x2 = min(outage(i,2), xl(2));
        if x2 <= x1, continue; end          % окно вне диапазона данных
        xv = [x1 x2 x2 x1];
        yv = [yl(1) yl(1) yl(2) yl(2)];
        p = patch(ax, xv, yv, [1 0.86 0.89], 'EdgeColor','none', ...
                  'HandleVisibility','off');
        uistack(p, 'bottom');               % под кривые
    end

    % Восстанавливаем пределы: патч не должен их менять
    xlim(ax, xl);
    ylim(ax, yl);
    set(ax, 'Layer', 'top');                % сетка поверх заливки
end