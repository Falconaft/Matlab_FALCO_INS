function plot_imu_comparison(results, cep_target)
%PLOT_IMU_COMPARISON  Сравнительные графики кандидатов IMU (Шаг 4).
%
%   Строит три фигуры:
%     (1) CDF горизонтальной ошибки по кандидатам (главный график сравнения)
%     (2) Scatter точек попадания в ENU с кругами КВО
%     (3) Столбчатая диаграмма КВО и R95 по кандидатам
%
%   Вход:
%     results     - cell-массив структур mc (из run_montecarlo_coast)
%     cep_target  - целевой КВО [м] для опорной линии

    n = numel(results);
    colors = lines(n);

    % =====================================================================
    % ФИГУРА 1: CDF горизонтальной ошибки (кумулятивная функция)
    % =====================================================================
    figure('Color','w','Position',[100 100 800 550]);
    hold on;
    for j = 1:n
        e = sort(results{j}.err_horiz);
        cdf = (1:numel(e))'/numel(e);            % эмпирическая CDF
        plot(e, cdf*100, 'LineWidth', 2.5, 'Color', colors(j,:), ...
             'DisplayName', results{j}.profile_name);
    end
    % Опорные линии: целевой КВО и уровень 50% (медиана = КВО)
    xline(cep_target, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Цель %d м', cep_target));
    yline(50, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
    grid on; box on;
    xlabel('Горизонтальная ошибка в точке удара [м]', 'FontSize', 11);
    ylabel('Доля прогонов ≤ ошибки [%]', 'FontSize', 11);
    title('CDF горизонтальной ошибки (чистая ИНС, без ESKF)', 'FontSize', 13, 'FontWeight','bold');
    legend('Location','southeast','FontSize',10);
    xlim([0, prctile(vertcat_field(results,'err_horiz'), 98)]);

    % =====================================================================
    % ФИГУРА 2: Scatter точек попадания в ENU с кругами КВО
    % =====================================================================
    figure('Color','w','Position',[920 100 700 700]);
    hold on;
    for j = 1:n
        E = results{j}.err_enu(:,1);
        Nn = results{j}.err_enu(:,2);
        scatter(E, Nn, 18, colors(j,:), 'filled', 'MarkerFaceAlpha', 0.4, ...
                'DisplayName', results{j}.profile_name);
        % Круг КВО (радиус = медиана горизонтальной ошибки)
        draw_circle(0, 0, results{j}.cep, colors(j,:));
    end
    % Целевой круг КВО
    draw_circle(0, 0, cep_target, [0 0 0]);
    plot(0,0,'k+','MarkerSize',14,'LineWidth',2,'HandleVisibility','off');
    grid on; box on; axis equal;
    xlabel('Ошибка Восток [м]', 'FontSize', 11);
    ylabel('Ошибка Север [м]', 'FontSize', 11);
    title('Разброс попаданий и круги КВО (чёрный = цель)', 'FontSize', 13, 'FontWeight','bold');
    legend('Location','best','FontSize',9);

    % =====================================================================
    % ФИГУРА 3: Столбчатая диаграмма КВО и R95
    % =====================================================================
    figure('Color','w','Position',[100 100 800 500]);
    cep_vals = cellfun(@(m) m.cep, results);
    r95_vals = cellfun(@(m) m.cep_r95, results);
    names = cellfun(@(m) m.profile_name, results, 'UniformOutput', false);
    bar_data = [cep_vals(:), r95_vals(:)];
    b = bar(bar_data, 'grouped');
    b(1).FaceColor = [0.2 0.5 0.8];
    b(2).FaceColor = [0.8 0.4 0.3];
    hold on;
    yline(cep_target, 'k--', 'LineWidth', 2, 'Label', sprintf('Цель КВО %d м', cep_target));
    set(gca, 'XTickLabel', names, 'FontSize', 10);
    ylabel('Ошибка [м]', 'FontSize', 11);
    title('КВО и 95%-й процентиль по кандидатам', 'FontSize', 13, 'FontWeight','bold');
    legend({'КВО (медиана)','R95'}, 'Location','best','FontSize',10);
    grid on;
end


% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% =========================================================================
function draw_circle(cx, cy, r, col)
%DRAW_CIRCLE  Рисует окружность радиуса r с центром (cx,cy).
    th = linspace(0, 2*pi, 100);
    plot(cx + r*cos(th), cy + r*sin(th), '--', 'Color', col, ...
         'LineWidth', 1.5, 'HandleVisibility','off');
end

function v = vertcat_field(results, field)
%VERTCAT_FIELD  Собирает поле field из всех результатов в один вектор.
    v = [];
    for j = 1:numel(results)
        v = [v; results{j}.(field)];  %#ok<AGROW>
    end
end
