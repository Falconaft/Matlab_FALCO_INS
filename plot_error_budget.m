function plot_error_budget(budget)
%PLOT_ERROR_BUDGET  Графики бюджета ошибок.
%
%   Две фигуры:
%     1. Столбчатая диаграмма вкладов + доля от общего КВО (парето)
%     2. КВО при отключении каждого источника (что даёт улучшение железа)

    n = numel(budget.names);
    contrib = budget.contrib;
    cepb = budget.cep_base;

    % =====================================================================
    % ФИГУРА 1: вклады источников (парето)
    % =====================================================================
    figure('Color','w','Position',[60 60 1000 560]);

    ax = subplot(1,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    barh(ax, contrib, 'FaceColor',[0.25 0.5 0.75], 'EdgeColor','w');
    set(ax, 'YTick', 1:n, 'YTickLabel', budget.names, 'YDir','reverse');
    xlabel(ax, 'Вклад в КВО, м');
    title(ax, sprintf('Вклад источников (КВО = %.2f м)', cepb), 'FontWeight','bold');
    % Подписи значений
    for i = 1:n
        text(ax, contrib(i)+0.02*max(contrib), i, sprintf('%.2f', contrib(i)), ...
             'VerticalAlignment','middle', 'FontSize',8);
    end

    ax = subplot(1,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    frac = 100*contrib.^2 / max(cepb^2, eps);
    cum  = cumsum(frac);
    bar(ax, frac, 'FaceColor',[0.35 0.65 0.45], 'EdgeColor','w');
    plot(ax, 1:n, cum, '-o', 'Color',[0.8 0.3 0.2], 'LineWidth',1.8, 'MarkerSize',5);
    yline(ax, 80, '--', 'Color',[0.5 0.5 0.5], 'Label','80%');
    set(ax, 'XTick', 1:n, 'XTickLabel', budget.names);
    xtickangle(ax, 40);
    ylabel(ax, 'Доля дисперсии, %');
    title(ax, 'Парето: накопленная доля','FontWeight','bold');
    legend(ax, {'вклад','накопленно'}, 'Location','east');

    % =====================================================================
    % ФИГУРА 2: КВО при отключении источника
    % =====================================================================
    figure('Color','w','Position',[100 100 900 500]);
    ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    barh(ax, budget.cep_off, 'FaceColor',[0.55 0.65 0.8], 'EdgeColor','w');
    xline(ax, cepb, 'r--', 'LineWidth',2, 'Label',sprintf('база %.2f м',cepb));
    set(ax, 'YTick', 1:n, 'YTickLabel', budget.names, 'YDir','reverse');
    xlabel(ax, 'КВО при отключённом источнике, м');
    title(ax, 'КВО при отключении каждого источника','FontWeight','bold');

    for i = 1:n
    text(ax, budget.cep_off(i)+0.02*cepb, i, ...
         sprintf('%.2f', budget.cep_off(i)), ...
         'VerticalAlignment','middle', 'FontSize',8);
    end
end
