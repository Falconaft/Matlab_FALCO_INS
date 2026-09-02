function m = scale_factor_metrics(s_true, s_est)
%SCALE_FACTOR_METRICS  Качество оценки масштабных коэффициентов (в ppm).
%
%   Тонкая обёртка над estimation_quality_metrics: та же математика, но
%   с переводом в ppm, поскольку именно так масштаб задан в datasheet.
%
%   КОНВЕНЦИЯ ОШИБКИ. Измерение искажается как v_изм = (1 + s_true)·v_ист,
%   механизация компенсирует оценкой: v_corr = (I - diag(s_est))·v_изм.
%   Остаточная относительная ошибка масштаба равна s_true - s_est, поэтому
%   ошибкой оценки считается err = s_est - s_true.
%
%   Вход:  s_true, s_est - (N x 3) безразмерные
%   Выход: m - метрики в ppm (см. estimation_quality_metrics)

    m = estimation_quality_metrics(s_true, s_est, 1e6);
end
