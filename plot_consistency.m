function plot_consistency(mc, p, cfg)
%PLOT_CONSISTENCY  Состоятельность фильтра: ANEES(t) и ANIS(t).
%
%   ЧТО ЭТО ПОКАЗЫВАЕТ. Метрики отвечают не на вопрос «насколько точен
%   фильтр», а на вопрос «честно ли он оценивает СОБСТВЕННУЮ
%   неопределённость». Фильтр может давать малую ошибку и при этом быть
%   несостоятельным — и наоборот.
%
%   ПОЧЕМУ УСРЕДНЕНИЕ ИДЁТ ПО АНСАМБЛЮ, А НЕ ПО ВРЕМЕНИ.
%   Отсчёты NEES/NIS внутри ОДНОЙ реализации КОРРЕЛИРОВАНЫ: состояние
%   фильтра меняется плавно, соседние моменты почти повторяют друг друга.
%   Поэтому среднее по времени нельзя трактовать как выборку независимых
%   chi-квадрат величин, и строить для него chi-квадрат коридор некорректно.
%
%   Классическая формулировка (Bar-Shalom) усредняет ПОПЕРЁК независимых
%   реализаций при ФИКСИРОВАННОМ времени:
%
%       ANEES(k) = (1/N) * sum_i NEES_i(t_k)
%       ANIS(k)  = (1/N) * sum_i NIS_i(t_k)
%
%   При согласованном фильтре sum_i ~ chi2(N*dim), поэтому коридор:
%
%       [ chi2(0.025, N*dim)/N , chi2(0.975, N*dim)/N ]
%
%   Идеальное значение — dim = 6 (позиция+скорость для ANEES, размерность
%   измерения GNSS для ANIS).
%
%   ИНТЕРПРЕТАЦИЯ:
%     выше коридора -> ПЕРЕУВЕРЕН: реальные ошибки больше ожидаемых,
%                      коридор +-3sigma занижен;
%     ниже коридора -> НЕДОУВЕРЕН: ковариация завышена, оценка неоптимальна.
%
%   Вход:
%     mc  - результат run_montecarlo_eskf2 (поля anees*, anis*)
%     p   - параметры траектории (опц., для отметки коаста)
%     cfg - конфигурация (опц., для отметки окна GNSS)

    if nargin < 2, p   = struct(); end
    if nargin < 3, cfg = struct(); end

    if ~isfield(mc,'anees') || isempty(mc.anees)
        warning('plot_consistency:noData', 'Метрики состоятельности отсутствуют.');
        return;
    end

    dim = mc.consist_dim;

    % Границы участков для подсветки
    t_boost_end = get_or(cfg, 'outage', [0 0; 0 0]);
    if isfield(cfg,'outage') && size(cfg.outage,1) >= 2
        t_gnss_on   = cfg.outage(1,2);
        t_coast_beg = cfg.outage(2,1);
    else
        t_gnss_on   = NaN;
        t_coast_beg = NaN;
    end

    figure('Color','w','Position',[100 100 1150 700]);

    % =====================================================================
    % (а) ANEES(t) — главный график
    % =====================================================================
    ax = subplot(2,2,[1 2]); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    t = mc.anees_t;
    plot(ax, t, mc.anees,    'LineWidth',1.5, 'Color',[0.15 0.35 0.70]);
    plot(ax, t, mc.anees_hi, '--', 'Color',[0.80 0.20 0.20], 'LineWidth',1.2);
    plot(ax, t, mc.anees_lo, '--', 'Color',[0.80 0.20 0.20], 'LineWidth',1.2);
    yline(ax, dim, 'k-', 'LineWidth',1.8);
    set(ax,'YScale','log');
    xlabel(ax,'Время, с'); ylabel(ax,'ANEES (лог. шкала)');
    title(ax, sprintf(['ANEES(t): усреднение по ансамблю N=%d при фиксированном ' ...
           'времени (идеал = %d)'], mc.consist_n_typ, dim), 'FontWeight','bold');
    legend(ax, {'ANEES','границы 95%','', sprintf('идеал = %d',dim)}, ...
           'Location','best');
    mark_phases(ax, t_gnss_on, t_coast_beg);

    % =====================================================================
    % (б) ANIS(t)
    % =====================================================================
    ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if isfield(mc,'anis') && ~isempty(mc.anis)
        tg = mc.anis_t;
        plot(ax, tg, mc.anis,    'LineWidth',1.4, 'Color',[0.20 0.55 0.35]);
        plot(ax, tg, mc.anis_hi, '--', 'Color',[0.80 0.20 0.20], 'LineWidth',1.1);
        plot(ax, tg, mc.anis_lo, '--', 'Color',[0.80 0.20 0.20], 'LineWidth',1.1);
        yline(ax, dim, 'k-', 'LineWidth',1.6);
        set(ax,'YScale','log');
        xlabel(ax,'Время, с'); ylabel(ax,'ANIS (лог. шкала)');
        title(ax,'ANIS(t) — вычислим В ПОЛЁТЕ без истины','FontWeight','bold');
        mark_phases(ax, t_gnss_on, t_coast_beg);
    end

    % =====================================================================
    % (в) Доля времени вне коридора
    % =====================================================================
    ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    above_e = mean(mc.anees > mc.anees_hi, 'omitnan')*100;
    below_e = mean(mc.anees < mc.anees_lo, 'omitnan')*100;
    inside_e = 100 - above_e - below_e;
    if isfield(mc,'anis') && ~isempty(mc.anis)
        above_n = mean(mc.anis > mc.anis_hi, 'omitnan')*100;
        below_n = mean(mc.anis < mc.anis_lo, 'omitnan')*100;
    else
        above_n = NaN; below_n = NaN;
    end
    inside_n = 100 - above_n - below_n;

    vals = [inside_e, above_e, below_e; inside_n, above_n, below_n];
    b = bar(ax, vals, 'stacked');
    b(1).FaceColor = [0.35 0.70 0.40];   % внутри
    b(2).FaceColor = [0.85 0.30 0.25];   % выше (переуверен)
    b(3).FaceColor = [0.35 0.45 0.80];   % ниже (недоуверен)
    set(ax,'XTick',1:2,'XTickLabel',{'ANEES','ANIS'});
    ylabel(ax,'Доля моментов, %'); ylim(ax,[0 100]);
    title(ax,'Где метрика лежит относительно коридора','FontWeight','bold');
    legend(ax, {'внутри 95%','выше (переуверен)','ниже (недоуверен)'}, ...
           'Location','best');
end


function mark_phases(ax, t_gnss_on, t_coast_beg)
%MARK_PHASES  Подсветка участков БЕЗ искажения пределов осей.
    xl = xlim(ax);   yl = ylim(ax);
    if ~isnan(t_gnss_on) && t_gnss_on > xl(1)
        pt = patch(ax, [xl(1) t_gnss_on t_gnss_on xl(1)], ...
                   [yl(1) yl(1) yl(2) yl(2)], [1 0.86 0.89], ...
                   'EdgeColor','none','HandleVisibility','off');
        uistack(pt,'bottom');
    end
    if ~isnan(t_coast_beg) && t_coast_beg < xl(2)
        pt = patch(ax, [t_coast_beg xl(2) xl(2) t_coast_beg], ...
                   [yl(1) yl(1) yl(2) yl(2)], [1 0.86 0.89], ...
                   'EdgeColor','none','HandleVisibility','off');
        uistack(pt,'bottom');
    end
    xlim(ax, xl);  ylim(ax, yl);
    set(ax,'Layer','top');
end


function v = get_or(s, name, default)
%GET_OR  Чтение поля структуры со значением по умолчанию.
    if isfield(s, name), v = s.(name); else, v = default; end
end
