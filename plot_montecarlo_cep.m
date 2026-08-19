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
    % ФИГУРА 4: качество оценки bias акселерометра
    % =====================================================================
    figure('Color','w','Position',[180 180 950 430]);

    ax = subplot(1,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    scatter(ax, mc.ba_true/g0*1e3, mc.ba_resid/g0*1e3, 26, [0.2 0.6 0.35], 'filled', ...
            'MarkerFaceAlpha',0.5);
    xlabel(ax,'Истинный turn-on bias, мг'); ylabel(ax,'Остаточная ошибка оценки, мг');
    title(ax,'Качество оценки bias акселерометра','FontWeight','bold');

    ax = subplot(1,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    frac = 100*(1 - mc.ba_resid./mc.ba_true);
    histogram(ax, frac, max(10,round(sqrt(numel(frac)))), 'FaceColor',[0.2 0.6 0.35], ...
              'EdgeColor','w');
    xline(ax, median(frac), 'r--', 'LineWidth',2, 'Label',sprintf('медиана %.1f%%',median(frac)));
    xlabel(ax,'Доля оценённого bias, %'); ylabel(ax,'Число прогонов');
    title(ax,'Насколько фильтр "выучил" bias','FontWeight','bold');
end


function draw_circle(ax, r, col, lw)
%DRAW_CIRCLE  Окружность радиуса r с центром в нуле.
    th = linspace(0, 2*pi, 200);
    plot(ax, r*cos(th), r*sin(th), '--', 'Color', col, 'LineWidth', lw);
end
