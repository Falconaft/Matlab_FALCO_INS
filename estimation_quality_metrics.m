function m = estimation_quality_metrics(v_true, v_est, unit_k)
%ESTIMATION_QUALITY_METRICS  Качество оценки параметра фильтром.
%
%   Общая функция для любых оцениваемых параметров (bias, масштаб и т.д.).
%   Введена, чтобы формулы жили в одном месте: scale_factor_metrics —
%   тонкая обёртка над ней.
%
%   ОСНОВНЫЕ ПОКАЗАТЕЛИ — RMS истинного значения и RMSE остатка. Именно они
%   информативны сами по себе: RMS показывает масштаб задачи, RMSE — что
%   осталось после работы фильтра.
%
%   ВСПОМОГАТЕЛЬНЫЙ показатель устранённой доли:
%       eta = 100 * (1 - RMSE / RMS)
%   Он удобен для сравнения, но НЕУСТОЙЧИВ при малом RMS: деление на малое
%   число раздувает результат. Поэтому eta следует читать только вместе
%   с RMS и RMSE, а не вместо них.
%
%   Интерпретация eta:
%       ~ 100%   оценка практически полная;
%       ~   0%   бесполезна (эквивалентна нулевой оценке);
%       <   0    ВРЕДНА: внесено больше ошибки, чем убрано.
%
%   Вход:
%     v_true - (N x 3) истинные значения
%     v_est  - (N x 3) оценки фильтра (те же единицы)
%     unit_k - множитель перевода в целевые единицы (например 1e6 для ppm)
%   Выход:
%     m - структура: .rms_true_axis, .rmse_axis, .eta_axis (1x3)
%                    .rms_true, .rmse, .eta (сводно)
%                    .n — число учтённых реализаций

    if nargin < 3 || isempty(unit_k), unit_k = 1; end

    ok = all(isfinite(v_true), 2) & all(isfinite(v_est), 2);
    vt = v_true(ok, :) * unit_k;
    ve = v_est(ok, :)  * unit_k;
    m.n = sum(ok);

    if m.n == 0
        m.rms_true_axis = nan(1,3);  m.rmse_axis = nan(1,3);
        m.eta_axis = nan(1,3);
        m.rms_true = NaN;  m.rmse = NaN;  m.eta = NaN;
        return;
    end

    err = ve - vt;

    m.rms_true_axis = sqrt(mean(vt.^2,  1));
    m.rmse_axis     = sqrt(mean(err.^2, 1));
    m.eta_axis      = 100 * (1 - m.rmse_axis ./ max(m.rms_true_axis, eps));

    m.rms_true = sqrt(mean(vt(:).^2));
    m.rmse     = sqrt(mean(err(:).^2));
    m.eta      = 100 * (1 - m.rmse / max(m.rms_true, eps));
end
