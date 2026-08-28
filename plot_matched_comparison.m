function plot_matched_comparison(mc_base, mc_match, cfg)
%PLOT_MATCHED_COMPARISON  Baseline против согласованной модели.
%
%   Сравнивает состоятельность фильтра в двух условиях при НЕИЗМЕННОМ
%   фильтре: в baseline истина содержит источники, отсутствующие в модели
%   фильтра; в matched-model они из истины убраны.
%
%   Если несостоятельность вызвана именно этими источниками, в matched-model
%   ANEES должен приблизиться к идеалу. Если нет — причина в другом
%   (например, в самой настройке P/Q).
%
%   Панели:
%     ANEES6(t), ANEES_r(t), ANEES_v(t), ANIS(t) — обе кампании вместе
%     с коридором 95% и линией идеала;
%     отношение diag(H·P·H')/diag(R) — проверка доминирования R.
%
%   Вход: mc_base, mc_match - результаты run_montecarlo_eskf2
%         cfg               - конфигурация (для отметки участков)

    if nargin < 3, cfg = struct(); end

    if isfield(cfg,'outage') && size(cfg.outage,1) >= 2
        t_gnss_on   = cfg.outage(1,2);
        t_coast_beg = cfg.outage(2,1);
    else
        t_gnss_on = NaN;  t_coast_beg = NaN;
    end

    cB = [0.80 0.25 0.20];      % baseline
    cM = [0.15 0.45 0.75];      % matched

    figure('Color','w','Position',[80 80 1200 760]);

    % --- (а) ANEES6 ---
    ax = subplot(3,2,[1 2]); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot_pair(ax, mc_base, mc_match, 'anees', 'anees_hi', 'anees_lo', ...
              6, cB, cM, t_gnss_on, t_coast_beg);
    ylabel(ax,'ANEES6'); title(ax,'ANEES6(t): позиция + скорость (идеал 6)', ...
           'FontWeight','bold');
    legend(ax, {'baseline','matched','коридор 95%','','идеал'}, 'Location','best');

    % --- (б) ANEES_r ---
    ax = subplot(3,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(mc_base,'anees_r')
        plot_pair(ax, mc_base, mc_match, 'anees_r', 'anees_r_hi', 'anees_r_lo', ...
                  3, cB, cM, t_gnss_on, t_coast_beg);
    end
    ylabel(ax,'ANEES_r'); title(ax,'Канал ПОЗИЦИИ (идеал 3)','FontWeight','bold');

    % --- (в) ANEES_v ---
    ax = subplot(3,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(mc_base,'anees_v')
        plot_pair(ax, mc_base, mc_match, 'anees_v', 'anees_v_hi', 'anees_v_lo', ...
                  3, cB, cM, t_gnss_on, t_coast_beg);
    end
    ylabel(ax,'ANEES_v'); title(ax,'Канал СКОРОСТИ (идеал 3)','FontWeight','bold');

    % --- (г) ANIS ---
    ax = subplot(3,2,5); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot_pair(ax, mc_base, mc_match, 'anis', 'anis_hi', 'anis_lo', ...
              6, cB, cM, t_gnss_on, t_coast_beg, 'anis_t');
    xlabel(ax,'Время, с'); ylabel(ax,'ANIS');
    title(ax,'ANIS(t) — измерение GNSS (идеал 6)','FontWeight','bold');

    % --- (д) вклад P против R ---
    ax = subplot(3,2,6); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(mc_base,'rho_pos')
        plot(ax, mc_base.rho_t,  mc_base.rho_pos,  '-',  'LineWidth',1.4, 'Color',cB);
        plot(ax, mc_base.rho_t,  mc_base.rho_vel,  '--', 'LineWidth',1.2, 'Color',cB);
        plot(ax, mc_match.rho_t, mc_match.rho_pos, '-',  'LineWidth',1.4, 'Color',cM);
        plot(ax, mc_match.rho_t, mc_match.rho_vel, '--', 'LineWidth',1.2, 'Color',cM);
        yline(ax, 1, 'k--', 'LineWidth',1.4, 'Label','P = R');
        set(ax,'YScale','log');
        xlabel(ax,'Время, с'); ylabel('\rho = diag(HPH'')/diag(R)');
        title(ax,'Доминирует ли R (сплошная — позиция, штрих — скорость)', ...
              'FontWeight','bold');
        legend(ax, {'base поз.','base скор.','match поз.','match скор.'}, ...
               'Location','best');
    end
end


function plot_pair(ax, mB, mM, fld, fhi, flo, ideal, cB, cM, tg, tc, tfld)
%PLOT_PAIR  Две кампании на одной панели с коридором и идеалом.
    if nargin < 12, tfld = 'anees_t'; end
    if ~isfield(mB, fld), return; end

    tB = mB.(tfld);   tM = mM.(tfld);
    plot(ax, tB, mB.(fld), 'LineWidth',1.5, 'Color',cB);
    plot(ax, tM, mM.(fld), 'LineWidth',1.5, 'Color',cM);
    plot(ax, tB, mB.(fhi), '--', 'Color',[0.45 0.45 0.45], 'LineWidth',1.0);
    plot(ax, tB, mB.(flo), '--', 'Color',[0.45 0.45 0.45], 'LineWidth',1.0);
    yline(ax, ideal, 'k-', 'LineWidth',1.6);
    set(ax,'YScale','log');

    % Подсветка участков без GNSS — после данных, с восстановлением пределов
    xl = xlim(ax);  yl = ylim(ax);
    for seg = {[xl(1) tg], [tc xl(2)]}
        s = seg{1};
        if all(isfinite(s)) && s(2) > s(1)
            pt = patch(ax, [s(1) s(2) s(2) s(1)], [yl(1) yl(1) yl(2) yl(2)], ...
                       [1 0.86 0.89], 'EdgeColor','none','HandleVisibility','off');
            uistack(pt,'bottom');
        end
    end
    xlim(ax, xl);  ylim(ax, yl);
    set(ax,'Layer','top');
end
