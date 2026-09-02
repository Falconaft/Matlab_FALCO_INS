function plot_scale_diagnostics(res, truth, imu, e_imu, p, cfg)
%PLOT_SCALE_DIAGNOSTICS  Диагностика оценки масштабных коэффициентов.
%
%   ЗАЧЕМ. Масштабные коэффициенты входят в состояние как СЛУЧАЙНЫЕ
%   КОНСТАНТЫ (F = 0, Qc = 0), то есть их ковариация может только убывать.
%   Наблюдаемость при этом целиком определяется ВОЗБУЖДЕНИЕМ: ошибка от
%   scale пропорциональна самому сигналу,
%       δa = -C·diag(f^b)·δs_a,      δψ̇ = -C·diag(ω^b)·δs_g
%   поэтому ось, вдоль которой сигнал мал, наблюдаться не может в принципе.
%   Панели возбуждения показаны рядом именно для этого: они объясняют,
%   почему конкретная ось оценивается хорошо или плохо.
%
%   Две фигуры (гироскоп и акселерометр), в каждой по три оси:
%     сплошная  — накопленная оценка фильтра после feedback/injection;
%     пунктир   — истинное значение;
%     заливка   — коридор ±3σ вокруг ИСТИННОГО значения (из diag(P));
%     снизу     — возбуждающий сигнал по той же оси.
%
%   Все масштабные величины в ppm.
%
%   Вход:
%     res    - результат eskf_run (нужны sg_est, sa_est, sig_sg, sig_sa)
%     truth  - истинная траектория (для сетки времени)
%     imu    - ИДЕАЛЬНЫЕ показания (.fb, .wib_b) — сигнал возбуждения
%     e_imu  - реализация ошибок (.gyro_scale, .accel_scale) — истина
%     p      - параметры траектории (опц., для отметки коаста)
%     cfg    - конфигурация (опц., для отметки окна GNSS)

    if nargin < 5, p   = struct(); end
    if nargin < 6, cfg = struct(); end

    if ~isfield(res,'sg_est')
        warning('plot_scale_diagnostics:noData', ...
                'В res нет полей sg_est/sa_est. Перезапусти eskf_run.');
        return;
    end

    PPM = 1e6;
    t   = res.t;
    ax_lbl = {'X','Y','Z'};

    % Границы участков
    if isfield(cfg,'outage') && size(cfg.outage,1) >= 2
        t_gnss_on   = cfg.outage(1,2);
        t_coast_beg = cfg.outage(2,1);
    else
        t_gnss_on = NaN;  t_coast_beg = NaN;
    end

    % Возбуждение прореживаем до сетки диагностики, чтобы графики были лёгкими
    exc_g = interp1(imu.t, imu.wib_b, t, 'linear');      % рад/с
    exc_a = interp1(imu.t, imu.fb,    t, 'linear');      % м/с²

    % --- Гироскоп ---
    draw_sensor('Гироскоп', 'sg_est', 'sig_sg', e_imu.gyro_scale, ...
                exc_g, '|\omega^b|, °/с', 180/pi);

    % --- Акселерометр ---
    draw_sensor('Акселерометр', 'sa_est', 'sig_sa', e_imu.accel_scale, ...
                exc_a, '|f^b|, g', 1/9.80665);


    % =====================================================================
    function draw_sensor(name, fld_est, fld_sig, s_true, exc, exc_lbl, exc_k)
        figure('Color','w','Position',[70 70 1150 700]);

        s_est = res.(fld_est) * PPM;        % ppm
        s_sig = res.(fld_sig) * PPM;        % ppm
        s_tru = s_true(:)'    * PPM;        % 1x3, ppm

        for k = 1:3
            % ---- Оценка и коридор ----
            axk = subplot(3,2,2*k-1); hold(axk,'on'); grid(axk,'on'); box(axk,'on');

            % Коридор ±3σ строится ВОКРУГ ИСТИННОГО значения: так сразу
            % видно, накрывает ли собственная неопределённость фильтра
            % фактическое значение параметра.
            hi = s_tru(k) + 3*s_sig(:,k);
            lo = s_tru(k) - 3*s_sig(:,k);
            fill(axk, [t; flipud(t)], [hi; flipud(lo)], [0.85 0.90 0.97], ...
                 'EdgeColor','none', 'HandleVisibility','off');

            plot(axk, t, s_est(:,k), 'LineWidth',1.6, 'Color',[0.15 0.35 0.70]);
            yline(axk, s_tru(k), 'r--', 'LineWidth',1.5);

            ylabel(axk, sprintf('s_%s, ppm', ax_lbl{k}));
            if k == 1
                title(axk, sprintf(['%s: оценка масштаба (пунктир — истина, ' ...
                       'заливка — \\pm3\\sigma)'], name), 'FontWeight','bold');
                legend(axk, {'оценка','истина'}, 'Location','best');
            end
            if k == 3, xlabel(axk,'Время, с'); end
            mark_phases(axk, t_gnss_on, t_coast_beg);

            % ---- Ошибка оценки ----
            axk = subplot(3,2,2*k); hold(axk,'on'); grid(axk,'on'); box(axk,'on');
            err = s_tru(k) - s_est(:,k);        % s_true - s_est
            plot(axk, t, err, 'LineWidth',1.5, 'Color',[0.75 0.30 0.20]);
            plot(axk, t,  3*s_sig(:,k), '--', 'Color',[0.45 0.45 0.45]);
            plot(axk, t, -3*s_sig(:,k), '--', 'Color',[0.45 0.45 0.45]);
            yline(axk, 0, 'k:', 'LineWidth',1.0);
            ylabel(axk, sprintf('ошибка s_%s, ppm', ax_lbl{k}));
            if k == 1
                title(axk, sprintf('Ошибка s_{true} - s_{est} и \\pm3\\sigma'), ...
                      'FontWeight','bold');
            end
            if k == 3, xlabel(axk,'Время, с'); end
            mark_phases(axk, t_gnss_on, t_coast_beg);
        end

        % ---- Возбуждение отдельной фигурой ----
        figure('Color','w','Position',[110 110 950 420]);
        axe = axes; hold(axe,'on'); grid(axe,'on'); box(axe,'on');
        for k = 1:3
            plot(axe, t, abs(exc(:,k))*exc_k, 'LineWidth',1.3);
        end
        set(axe,'YScale','log');
        xlabel(axe,'Время, с'); ylabel(axe, exc_lbl);
        title(axe, sprintf(['%s: ВОЗБУЖДЕНИЕ по осям — определяет ' ...
               'наблюдаемость масштаба'], name), 'FontWeight','bold');
        legend(axe, ax_lbl, 'Location','best');
        mark_phases(axe, t_gnss_on, t_coast_beg);
    end
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
