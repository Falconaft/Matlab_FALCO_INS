function plot_montecarlo_cep(mc, cep_target)
%PLOT_MONTECARLO_CEP  Графики результатов Монте-Карло кампании ESKF.
%
%   Четыре фигуры:
%     1. Разброс точек попадания в ENU с кругами КВО и R95
%     2. Гистограмма + CDF горизонтальной ошибки
%     3. Сходимость оценки КВО по мере роста выборки (сколько прогонов хватает)
%     4. Качество оценки bias акселерометра фильтром
%
%   Вход:
%     mc         - результат run_montecarlo_eskf
%     cep_target - целевой КВО [м] для опорной линии

    if nargin < 2, cep_target = 10; end

    ok = ~isnan(mc.err_horiz);
    h  = mc.err_horiz(ok);
    E  = mc.err_enu(ok,1);
    Nn = mc.err_enu(ok,2);
    g0 = 9.80665;

    % =====================================================================
    % ФИГУРА 1: разброс попаданий с кругами КВО
    % =====================================================================
    figure('Color','w','Position',[60 60 760 720]);
    ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');

    scatter(ax, E, Nn, 26, [0.15 0.4 0.75], 'filled', 'MarkerFaceAlpha',0.45);
    draw_circle(ax, mc.cep,     [0.85 0.2 0.2], 2.0);   % круг КВО
    draw_circle(ax, mc.r95,     [0.9 0.55 0.1], 1.5);   % круг R95
    draw_circle(ax, cep_target, [0 0 0],        2.0);   % целевой круг
    plot(ax, 0, 0, 'k+', 'MarkerSize',16, 'LineWidth',2);

    xlabel(ax,'Ошибка Восток, м'); ylabel(ax,'Ошибка Север, м');
    title(ax, sprintf('%s: разброс попаданий, N=%d\nКВО=%.2f м, R95=%.2f м (цель %.0f м)', ...
          mc.prof_name, numel(h), mc.cep, mc.r95, cep_target), 'FontWeight','bold');
    legend(ax, {'попадания', sprintf('КВО %.2f м',mc.cep), ...
                sprintf('R95 %.2f м',mc.r95), sprintf('цель %.0f м',cep_target)}, ...
           'Location','best');

    % =====================================================================
    % ФИГУРА 2: гистограмма и CDF
    % =====================================================================
    figure('Color','w','Position',[100 100 1000 440]);

    ax = subplot(1,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    histogram(ax, h, max(10,round(sqrt(numel(h)))), 'FaceColor',[0.3 0.55 0.8], ...
              'EdgeColor','w');
    xline(ax, mc.cep, 'r--', 'LineWidth',2, 'Label','КВО');
    xline(ax, mc.r95, 'Color',[0.9 0.55 0.1], 'LineStyle','--', 'LineWidth',1.5, 'Label','R95');
    xlabel(ax,'Горизонтальная ошибка, м'); ylabel(ax,'Число прогонов');
    title(ax,'Распределение ошибки','FontWeight','bold');

    ax = subplot(1,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    hs  = sort(h);
    cdf = (1:numel(hs))'/numel(hs)*100;
    plot(ax, hs, cdf, 'LineWidth',2.2, 'Color',[0.15 0.4 0.75]);
    yline(ax, 50, 'r--', 'LineWidth',1.4, 'Label','50% = КВО');
    yline(ax, 95, 'Color',[0.9 0.55 0.1], 'LineStyle','--', 'LineWidth',1.2, 'Label','95%');
    xline(ax, cep_target, 'k--', 'LineWidth',1.8);
    xlabel(ax,'Горизонтальная ошибка, м'); ylabel(ax,'Доля прогонов \leq ошибки, %');
    title(ax, sprintf('CDF (в цель уложилось %.1f%%)', 100*mc.frac_pass),'FontWeight','bold');
    ylim(ax,[0 100]);

    % =====================================================================
    % ФИГУРА 3: сходимость оценки КВО по размеру выборки
    % =====================================================================
    figure('Color','w','Position',[140 140 900 430]);
    ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    n_run  = numel(h);
    cep_run = zeros(n_run,1);
    for i = 1:n_run
        cep_run(i) = median(h(1:i));
    end
    plot(ax, 1:n_run, cep_run, 'LineWidth',1.8, 'Color',[0.15 0.4 0.75]);
    yline(ax, mc.cep, 'r--', 'LineWidth',1.5, 'Label','итоговый КВО');

    % Полоса неопределённости медианы ~ 1.253*sigma/sqrt(n)
    nn = (1:n_run)';
    band = 1.253*mc.std_horiz./sqrt(nn);
    plot(ax, nn, mc.cep+band, ':', 'Color',[0.5 0.5 0.5]);
    plot(ax, nn, mc.cep-band, ':', 'Color',[0.5 0.5 0.5]);

    xlabel(ax,'Число реализаций'); ylabel(ax,'Текущая оценка КВО, м');
    title(ax,'Сходимость оценки КВО (пунктир — ожидаемая точность медианы)','FontWeight','bold');

    % =====================================================================
    % ФИГУРА 4: качество оценки bias ПЕРЕД ТЕРМИНАЛЬНЫМ КОАСТОМ
    % =====================================================================
    % По осям: X — модуль ИСТИННОГО ПОЛНОГО смещения (turn-on + in-run + VRC),
    %          Y — модуль остатка оценки |b_true - b_est|.
    % Точки НИЖЕ диагонали означают, что фильтр уменьшил ошибку смещения;
    % выше — что оценка хуже, чем просто принять bias нулевым.
    %
    % Процент "сколько выучено" сознательно НЕ строится: при малом истинном
    % смещении отношение неустойчиво и даёт бессмысленные значения.
    figure('Color','w','Position',[180 180 1000 460]);
    deg = pi/180;

    % --- Акселерометр ---
    ax = subplot(1,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(mc,'ba_true_pre')
        xa = vecnorm(mc.ba_true_pre, 2, 2)/g0*1e3;      % мг
        ya = mc.ba_resid_pre/g0*1e3;                    % мг
        good = isfinite(xa) & isfinite(ya);
        scatter(ax, xa(good), ya(good), 30, [0.2 0.6 0.35], 'filled', ...
                'MarkerFaceAlpha',0.55);
        lim = [0, max([xa(good); ya(good); eps])*1.05];
        plot(ax, lim, lim, 'k--', 'LineWidth',1.2);
        xlim(ax, lim); ylim(ax, lim);
        xlabel(ax,'|истинное полное смещение|, мг');
        ylabel(ax,'|остаток оценки|, мг');
        title(ax, sprintf('Акселерометр (медианы %.3f / %.3f мг)', ...
              median(xa(good)), median(ya(good))), 'FontWeight','bold');
        legend(ax, {'реализации','y = x (оценка бесполезна)'}, 'Location','best');
    end

    % --- Гироскоп ---
    ax = subplot(1,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(mc,'bg_true_pre')
        xg = vecnorm(mc.bg_true_pre, 2, 2)/deg*3600;    % °/ч
        yg = mc.bg_resid_pre/deg*3600;                  % °/ч
        good = isfinite(xg) & isfinite(yg);
        scatter(ax, xg(good), yg(good), 30, [0.75 0.35 0.15], 'filled', ...
                'MarkerFaceAlpha',0.55);
        lim = [0, max([xg(good); yg(good); eps])*1.05];
        plot(ax, lim, lim, 'k--', 'LineWidth',1.2);
        xlim(ax, lim); ylim(ax, lim);
        xlabel(ax,'|истинное полное смещение|, °/ч');
        ylabel(ax,'|остаток оценки|, °/ч');
        title(ax, sprintf('Гироскоп (медианы %.3f / %.3f °/ч)', ...
              median(xg(good)), median(yg(good))), 'FontWeight','bold');
        legend(ax, {'реализации','y = x (оценка бесполезна)'}, 'Location','best');
    end
end


function draw_circle(ax, r, col, lw)
%DRAW_CIRCLE  Окружность радиуса r с центром в нуле.
    th = linspace(0, 2*pi, 200);
    plot(ax, r*cos(th), r*sin(th), '--', 'Color', col, 'LineWidth', lw);
end
