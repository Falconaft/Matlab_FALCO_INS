function mc = run_montecarlo_coast(profile, imu_ideal, truth, c, fnav, t_coast_start, N_mc, base_seed)
%RUN_MONTECARLO_COAST  Монте-Карло чистой ИНС на бесспутниковом коасте.
%
%   Для заданного профиля IMU прогоняет N_mc независимых реализаций:
%   на каждой - разыгрывает ошибки датчика, накладывает их на идеальные
%   показания, прогоняет ИНС-механизацию на 60-с коасте (старт из ИСТИННОГО
%   состояния, т.е. как если бы GNSS идеально поджал навигатор перед пропаданием),
%   и измеряет горизонтальную (в ENU) и полную радиальную ошибку в точке удара.
%
%   ВНИМАНИЕ: это оценка БЕЗ ESKF. Метрика показывает чистый инерциальный дрейф,
%   определяемый в основном turn-on bias. Финальный КВО (с оценкой bias фильтром)
%   считается после Шага 5. Здесь - валидация модели и первичный отсев.
%
%   Вход:
%     profile        - профиль IMU (из imu_profile_*.m)
%     imu_ideal      - идеальные показания (inverse_mech)
%     truth          - истинная траектория
%     c              - константы
%     fnav           - частота навигатора [Гц]
%     t_coast_start  - время начала коаста [с]
%     N_mc           - число прогонов Монте-Карло
%     base_seed      - базовый seed (для воспроизводимости)
%   Выход:
%     mc - структура с полями:
%          err_radial   (N_mc x 1) полная радиальная ошибка в точке удара [м]
%          err_horiz    (N_mc x 1) горизонтальная (ENU E-N) ошибка [м]
%          err_enu      (N_mc x 3) вектор ошибки в ENU [м]
%          cep          КВО (медиана горизонтальной ошибки) [м]
%          cep_r95      95-й процентиль горизонтальной ошибки [м]
%          profile_name имя профиля

    a0 = round(t_coast_start/truth.dt) + 1;
    r0 = truth.r0;
    Cen = truth.Cen;

    err_radial = zeros(N_mc,1);
    err_horiz  = zeros(N_mc,1);
    err_enu    = zeros(N_mc,3);

    for i = 1:N_mc
        % --- Независимый поток случайных чисел на каждый прогон ---
        rng_stream = RandStream('mt19937ar', 'Seed', base_seed + i);

        % --- Розыгрыш постоянных ошибок (turn-on bias, scale) ---
        e = imu_draw_errors(profile, rng_stream);

        % --- Наложение полной модели ошибок на идеальные показания ---
        imu_noisy = add_imu_errors(imu_ideal, profile, e, rng_stream);

        % --- Прогон навигатора на коасте от истинного состояния ---
        nav = ins_mechanize(imu_noisy, truth, c, fnav, a0);

        % --- Ошибка в точке удара (последняя точка навигатора) ---
        r_nav_end = nav.R(end,:)';
        r_true_end = truth.R(end,:)';
        d_ecef = r_nav_end - r_true_end;

        % Проекция ошибки в локальную ENU (для разделения горизонт/вертикаль)
        d_enu = Cen' * d_ecef;
        err_enu(i,:)   = d_enu';
        err_horiz(i)   = sqrt(d_enu(1)^2 + d_enu(2)^2);   % E-N горизонталь
        err_radial(i)  = norm(d_ecef);                    % полная 3D
    end

    % --- Статистика ---
    mc.err_radial   = err_radial;
    mc.err_horiz    = err_horiz;
    mc.err_enu      = err_enu;
    mc.cep          = median(err_horiz);                  % КВО = медиана горизонт. ошибки
    mc.cep_r95      = prctile(err_horiz, 95);
    mc.mean_horiz   = mean(err_horiz);
    mc.profile_name = profile.name;
end
