function plot_gyro_bias_diagnostics(res, imu, p, cfg)
%PLOT_GYRO_BIAS_DIAGNOSTICS  Диагностика оценки смещения гироскопа.
%
%   ЗАЧЕМ. В Шаге 9 медиана истинного полного смещения перед коастом
%   составляет ~1.0 °/ч, а остаток оценки ~5.4 °/ч — то есть оценка вносит
%   больше ошибки, чем убирает. Эти графики отвечают, какая ось портит
%   итоговую норму и почему.
%
%   ЧТО ПОКАЗЫВАЕТСЯ, по каждой оси:
%     истинное ПОЛНОЕ смещение (turn-on + in-run + VRC) — оно меняется
%       во времени, в отличие от одного лишь turn-on;
%     накопленная оценка фильтра bg_est после feedback/injection;
%     ошибка bg_true - bg_est с коридором ±3σ из diag(P(10:12,10:12)),
%       что сразу показывает, есть ли ПЕРЕУВЕРЕННОСТЬ по σ_bg;
%     корреляция δbg_i <-> δsg_i из самой P.
%
%   ПРО КОРРЕЛЯЦИЮ. Измерение гироскопа устроено так:
%       δω = δbg + diag(ω)·δsg
%   Оба состояния входят в одно и то же измерение, и разделить их можно
%   ТОЛЬКО за счёт изменения ω во времени: при постоянной ω они образуют
%   вырожденную комбинацию. Поэтому |rho| близкое к единице означает
%   неразделимость, а не ошибку фильтра.
%
%   Все величины смещения в °/ч.
%
%   Вход:
%     res - результат eskf_run (нужны bg_true, bg_est, sig_bg, rho_bg_sg)
%     imu - ИДЕАЛЬНЫЕ показания (.t, .wib_b) — возбуждение
%     p   - параметры траектории (опц.)
%     cfg - конфигурация (опц., для отметки участков)

    if nargin < 3, p   = struct(); end
    if nargin < 4, cfg = struct(); end

    if ~isfield(res,'bg_true')
        warning('plot_gyro_bias_diagnostics:noData', ...
                'В res нет поля bg_true. Перезапусти eskf_run.');
        return;
    end

    DPH = 180/pi*3600;          % рад/с -> °/ч
    DPS = 180/pi;               % рад/с -> °/с
    t   = res.t;
    ax_lbl = {'X','Y','Z'};

    if isfield(cfg,'outage') && size(cfg.outage,1) >= 2
        t_gnss_on   = cfg.outage(1,2);
        t_coast_beg = cfg.outage(2,1);
    else
        t_gnss_on = NaN;  t_coast_beg = NaN;
    end

    bg_tru = res.bg_true  * DPH;
    bg_est = res.bg_est   * DPH;
    bg_sig = res.sig_bg   * DPH;
    bg_err = bg_tru - bg_est;

    % =====================================================================
    % ФИГУРА 1: оценка и ошибка по осям
    % =====================================================================
    figure('Color','w','Position',[70 70 1180 720]);
    for k = 1:3
        % ---- Истина и оценка ----
        ax = subplot(3,2,2*k-1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        plot(ax, t, bg_tru(:,k), 'r--', 'LineWidth',1.5);
        plot(ax, t, bg_est(:,k), 'LineWidth',1.6, 'Color',[0.15 0.35 0.70]);
        ylabel(ax, sprintf('b_{g%s}, °/ч', ax_lbl{k}));
        if k == 1
            title(ax,'Смещение гироскопа: истина и оценка','FontWeight','bold');
            legend(ax, {'истина (полная)','оценка'}, 'Location','best');
        end
        if k == 3, xlabel(ax,'Время, с'); end
        mark_phases(ax, t_gnss_on, t_coast_beg);

        % ---- Ошибка и коридор ----
        ax = subplot(3,2,2*k); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        fill(ax, [t; flipud(t)], [3*bg_sig(:,k); flipud(-3*bg_sig(:,k))], ...
             [0.88 0.92 0.98], 'EdgeColor','none', 'HandleVisibility','off');
        plot(ax, t, bg_err(:,k), 'LineWidth',1.5, 'Color',[0.75 0.30 0.20]);
        yline(ax, 0, 'k:', 'LineWidth',1.0);
        ylabel(ax, sprintf('ошибка b_{g%s}, °/ч', ax_lbl{k}));
        if k == 1
            title(ax,'Ошибка b_{true} - b_{est}, заливка — \pm3\sigma', ...
                  'FontWeight','bold');
        end
        if k == 3, xlabel(ax,'Время, с'); end
        mark_phases(ax, t_gnss_on, t_coast_beg);
    end

    % =====================================================================
    % ФИГУРА 2: возбуждение и разделимость
    % =====================================================================
    figure('Color','w','Position',[110 110 1150 640]);

    % ---- Возбуждение по осям ----
    ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    w_exc = interp1(imu.t, imu.wib_b, t, 'linear') * DPS;
    for k = 1:3
        plot(ax, t, abs(w_exc(:,k)), 'LineWidth',1.3);
    end
    set(ax,'YScale','log');
    ylabel(ax,'|\omega^b|, °/с'); xlabel(ax,'Время, с');
    title(ax,'Возбуждение гироскопа по осям','FontWeight','bold');
    legend(ax, ax_lbl, 'Location','best');
    mark_phases(ax, t_gnss_on, t_coast_beg);

    % ---- Корреляция bg <-> sg ----
    ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(res,'rho_bg_sg')
        for k = 1:3
            plot(ax, t, res.rho_bg_sg(:,k), 'LineWidth',1.4);
        end
        yline(ax,  1, 'k--');  yline(ax, -1, 'k--');
        yline(ax,  0, 'k:');
        ylim(ax, [-1.05 1.05]);
        ylabel(ax,'\rho(\delta b_g, \delta s_g)'); xlabel(ax,'Время, с');
        title(ax,'Разделимость смещения и масштаба','FontWeight','bold');
        legend(ax, ax_lbl, 'Location','best');
        mark_phases(ax, t_gnss_on, t_coast_beg);
    end

    % ---- Ошибка против коридора: отношение ----
    ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    ratio = abs(bg_err) ./ max(bg_sig, eps);
    for k = 1:3
        plot(ax, t, ratio(:,k), 'LineWidth',1.3);
    end
    yline(ax, 3, 'r--', 'LineWidth',1.5, 'Label','3\sigma');
    set(ax,'YScale','log');
    ylabel(ax,'|ошибка| / \sigma'); xlabel(ax,'Время, с');
    title(ax,'Переуверенность: выход за 3\sigma','FontWeight','bold');
    legend(ax, ax_lbl, 'Location','best');
    mark_phases(ax, t_gnss_on, t_coast_beg);

    % ---- Сводка по осям в конце ----
    ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    k_end = numel(t);
    vals = [abs(bg_tru(k_end,:)); abs(bg_err(k_end,:)); 3*bg_sig(k_end,:)]';
    bar(ax, vals);
    set(ax,'XTick',1:3,'XTickLabel',ax_lbl);
    ylabel(ax,'°/ч');
    title(ax,'На конец сценария по осям','FontWeight','bold');
    legend(ax, {'|истина|','|ошибка|','3\sigma'}, 'Location','best');
end


function mark_phases(ax, t_gnss_on, t_coast_beg)
%MARK_PHASES  Подсветка участков без GNSS, без искажения пределов осей.
    xl = xlim(ax);   yl = ylim(ax);
    segs = {[xl(1) t_gnss_on], [t_coast_beg xl(2)]};
    for q = 1:2
        s = segs{q};
        if all(isfinite(s)) && s(2) > s(1)
            pt = patch(ax, [s(1) s(2) s(2) s(1)], [yl(1) yl(1) yl(2) yl(2)], ...
                       [1 0.86 0.89], 'EdgeColor','none', 'HandleVisibility','off');
            uistack(pt,'bottom');
        end
    end
    xlim(ax, xl);   ylim(ax, yl);
    set(ax,'Layer','top');
end
