function truth = gen_truth_v2(p, c)
%GEN_TRUTH_V2  Генератор реалистичной траектории (Шаг 7).
%
%   Динамика в ECEF с полной физикой:
%     - тяга через удельный импульс, масса убывает по расходу топлива
%     - стандартная атмосфера ISA, аэродинамика CD(M,alpha), CL(M,alpha)
%     - гравитация с J2, Кориолис, центробежное (4 члена механизации)
%     - фазовая программа угла атаки со СГЛАЖЕННЫМИ переходами
%
%   Ключевое отличие от gen_truth.m (Шаг 1): удельная сила НЕ обнуляется
%   после разгона. На баллистике |f| ~ 0.13..1.8 g вместо машинного нуля,
%   поэтому ошибка ориентации и scale factor становятся НАБЛЮДАЕМЫМИ.
%
%   Фазы:
%     1 - разгон (alpha = 0)
%     2 - перелом программы тангажа (alpha = p.alpha_ph2), пока gamma > theta_target
%     3 - после апогея, высота > h_switch (alpha = p.alpha_apogee)
%     4 - терминальный участок, высота < h_switch (alpha = p.alpha_low)
%
%   Вход:  p - параметры траектории (traj_params)
%          c - константы (constants)
%   Выход: truth - структура истинной траектории:
%            .t (N x 1) время [с]
%            .R (N x 3) позиция ECEF [м]
%            .V (N x 3) скорость ECEF [м/с]
%            .f_ecef (N x 3) удельная сила в ECEF [м/с²]
%            .alpha (N x 1) угол атаки [град]
%            .phase (N x 1) номер фазы
%            .alt, .Vmag, .gamma, .Mach, .mass, .thrust - диагностика
%            .dt, .r0, .Cen, .nb (индекс конца разгона)
%            .C_e_ned0 (3x3) NED -> ECEF в точке старта
%            .C_ned_e0 (3x3) ECEF -> NED в точке старта (= C_e_ned0')
%            .ned (N x 3) позиция в NED относительно точки старта [м]

    dt = p.dt;

    % =====================================================================
    % Шаг 7.1: НАЧАЛЬНЫЕ УСЛОВИЯ
    % =====================================================================
    r0  = lla2ecef(deg2rad(p.lat0), deg2rad(p.lon0), p.h0, c);
    Cen = C_e_n(deg2rad(p.lat0), deg2rad(p.lon0));   % столбцы: E, N, U

    % Матрица NED в точке старта. Конвенция та же, что у C_e_n:
    % столбцы — орты [Север, Восток, Вниз] в ECEF, применение v_e = C * v_ned.
    C_e_ned0 = C_e_ned(deg2rad(p.lat0), deg2rad(p.lon0));

    % Начальная скорость по азимуту и углу возвышения (в ENU, затем в ECEF)
    az = deg2rad(p.azimuth);
    el = deg2rad(p.launch_el);
    v_enu = [ p.V0*cos(el)*sin(az);      % Восток
              p.V0*cos(el)*cos(az);      % Север
              p.V0*sin(el) ];            % Верх
    v0 = Cen * v_enu;

    % =====================================================================
    % Шаг 7.2: ПРЕАЛЛОКАЦИЯ
    % =====================================================================
    N_max = ceil(p.t_max/dt) + 1;
    t      = zeros(N_max,1);
    R      = zeros(N_max,3);   V      = zeros(N_max,3);
    F_e    = zeros(N_max,3);   alt    = zeros(N_max,1);
    Vmag   = zeros(N_max,1);   gamma  = zeros(N_max,1);
    Mach   = zeros(N_max,1);   alpha  = zeros(N_max,1);
    mass   = zeros(N_max,1);   thrust = zeros(N_max,1);
    phase  = zeros(N_max,1);

    % =====================================================================
    % Шаг 7.3: ДИСКРЕТНОЕ СОСТОЯНИЕ ФАЗЫ
    % =====================================================================
    % ВАЖНО: фазовое состояние обновляется ВНЕ RK4 и постоянно на такте.
    % Иначе интегратор ловит разрывы между стадиями k1..k4.
    ds.phase        = 1;
    ds.ph2_done     = false;
    ds.passed_apogee= false;
    ds.max_alt      = -1e9;
    ds.alpha_cmd    = 0.0;      % целевой угол атаки текущей фазы
    ds.alpha_from   = 0.0;      % угол, от которого идёт перекладка
    ds.t_switch     = -1e9;     % момент начала перекладки

    % =====================================================================
    % Шаг 7.4: ОСНОВНОЙ ЦИКЛ ИНТЕГРИРОВАНИЯ (RK4)
    % =====================================================================
    y  = [r0; v0];
    tk = 0.0;
    k  = 0;

    while tk < p.t_max
        k = k + 1;

        % --- Оценка текущего состояния (для логов и фазовой логики) ---
        st = traj_state(tk, y, p, c, ds);

        % --- Обновление фазы (вне RK4!) ---
        ds = update_phase(tk, st, ds, p);

        % --- Запись в лог ---
        t(k)      = tk;
        R(k,:)    = y(1:3)';    V(k,:)   = y(4:6)';
        F_e(k,:)  = st.f_ecef';
        alt(k)    = st.alt;     Vmag(k)  = st.V;
        gamma(k)  = st.gamma;   Mach(k)  = st.Mach;
        alpha(k)  = st.alpha;   mass(k)  = st.mass;
        thrust(k) = st.thrust;  phase(k) = ds.phase;

        % --- Условие удара ---
        if st.alt <= p.impact_alt && tk > 10.0
            break;
        end

        % --- Шаг RK4 (фазовое состояние ds заморожено на такте) ---
        s1 = traj_state(tk,          y,               p, c, ds);
        s2 = traj_state(tk + dt/2,   y + dt/2*s1.dy,  p, c, ds);
        s3 = traj_state(tk + dt/2,   y + dt/2*s2.dy,  p, c, ds);
        s4 = traj_state(tk + dt,     y + dt*s3.dy,    p, c, ds);

        y  = y + dt/6*(s1.dy + 2*s2.dy + 2*s3.dy + s4.dy);
        tk = tk + dt;
    end

    N = k;

    % =====================================================================
    % Шаг 7.5: УПАКОВКА РЕЗУЛЬТАТА
    % =====================================================================
    truth.t      = t(1:N);
    truth.R      = R(1:N,:);       truth.V      = V(1:N,:);
    truth.f_ecef = F_e(1:N,:);
    truth.alt    = alt(1:N);       truth.Vmag   = Vmag(1:N);
    truth.gamma  = gamma(1:N);     truth.Mach   = Mach(1:N);
    truth.alpha  = alpha(1:N);     truth.mass   = mass(1:N);
    truth.thrust = thrust(1:N);    truth.phase  = phase(1:N);
    truth.dt     = dt;
    truth.r0     = r0;
    truth.Cen    = Cen;

    % --- NED-интерфейс (внутренний core остаётся ECEF) ---
    truth.C_e_ned0 = C_e_ned0;        % NED -> ECEF
    truth.C_ned_e0 = C_e_ned0';       % ECEF -> NED

    % Индекс конца работы двигателя (для диагностики)
    truth.nb = min(round(p.t_burn/dt) + 1, N);

    % Локальные ENU-координаты относительно точки старта
    truth.enu = (truth.R - r0') * Cen;

    % Позиция в NED относительно точки старта. Величины хранятся ПОСТРОЧНО,
    % поэтому преобразование — умножение справа на матрицу NED -> ECEF.
    % Связь с ENU: ned = [N E -U], то есть
    %   truth.ned(:,1) ==  truth.enu(:,2)
    %   truth.ned(:,2) ==  truth.enu(:,1)
    %   truth.ned(:,3) == -truth.enu(:,3)
    truth.ned = (truth.R - r0') * C_e_ned0;
end


% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% =========================================================================

function st = traj_state(tk, y, p, c, ds)
%TRAJ_STATE  Полное состояние в точке: силы, производные, диагностика.

    r = y(1:3);
    v = y(4:6);

    % --- Геодезические координаты и локальный базис ---
    [lat, lon, alt] = ecef2lla(r, c);
    Cen = C_e_n(lat, lon);            % столбцы: E, N, U

    % --- Скорость в локальной ENU и параметры траектории ---
    v_enu = Cen' * v;
    vE = v_enu(1);   vN = v_enu(2);   vU = v_enu(3);
    V  = max(0.1, norm(v_enu));
    v_horiz = hypot(vE, vN);
    gamma_rad = atan2(vU, max(1e-9, v_horiz));     % угол наклона траектории
    track_rad = atan2(vE, vN);                     % путевой угол

    % --- Атмосфера и число Маха ---
    [rho, a_sound] = atmosphere_isa(alt);
    Mach  = V / a_sound;
    q_inf = 0.5 * rho * V^2;                       % скоростной напор

    % --- Масса и тяга со СГЛАЖЕННОЙ отсечкой ---
    m_fuel = max(0.0, p.m_fuel - p.fuel_rate*min(tk, p.t_burn));
    m_curr = p.m_control + p.m_warhead + p.m_casing + m_fuel;

    if tk < p.t_burn
        thrust_frac = 1.0;
    else
        thrust_frac = 1.0 - smoothstep5( (tk - p.t_burn)/p.t_cutoff_smooth );
    end
    thrust_val = p.thrust * thrust_frac;

    % --- Угол атаки со СГЛАЖЕННОЙ перекладкой ---
    s_alpha   = smoothstep5( (tk - ds.t_switch)/p.t_alpha_smooth );
    alpha_deg = ds.alpha_from + (ds.alpha_cmd - ds.alpha_from)*s_alpha;
    alpha_rad = deg2rad(alpha_deg);

    % --- Аэродинамические силы ---
    [CD, CL] = aero_coeffs(Mach, alpha_rad);
    drag = CD * q_inf * p.area;
    lift = CL * q_inf * p.area;

    % --- Силы в скоростной системе (продольная и нормальная) ---
    F_tangent = thrust_val*cos(alpha_rad) - drag;
    F_normal  = thrust_val*sin(alpha_rad) + lift;

    % --- Удельная сила в ENU ---
    % Продольная - вдоль вектора скорости; нормальная - перпендикулярно ему
    % в вертикальной плоскости.
    cg = cos(gamma_rad);   sg = sin(gamma_rad);
    ct = cos(track_rad);   stk = sin(track_rad);

    f_N = ( F_tangent*cg*ct - F_normal*sg*ct ) / m_curr;
    f_E = ( F_tangent*cg*stk - F_normal*sg*stk ) / m_curr;
    f_U = ( F_tangent*sg + F_normal*cg ) / m_curr;

    f_enu  = [f_E; f_N; f_U];
    f_ecef = Cen * f_enu;

    % --- Полное ускорение: 4 члена механизации ECEF ---
    WIEx = skew(c.WIE);
    acc  = f_ecef + grav_ecef(r, c) - 2*WIEx*v - WIEx*(WIEx*r);

    % --- Упаковка ---
    st.dy     = [v; acc];
    st.f_ecef = f_ecef;
    st.alt    = alt;
    st.V      = V;
    st.gamma  = rad2deg(gamma_rad);
    st.Mach   = Mach;
    st.alpha  = alpha_deg;
    st.mass   = m_curr;
    st.thrust = thrust_val;
end


function ds = update_phase(tk, st, ds, p)
%UPDATE_PHASE  Обновление дискретного фазового состояния (вне RK4).

    % --- Детектор апогея ---
    if st.alt > ds.max_alt
        ds.max_alt = st.alt;
    elseif st.alt < ds.max_alt - 1.0
        ds.passed_apogee = true;
    end

    % --- Определение текущей фазы ---
    if ds.passed_apogee
        if st.alt < p.h_switch
            ds.phase = 4;
        else
            ds.phase = 3;
        end
    elseif tk >= p.t_ph2
        ds.phase = 2;
        % Перелом тангажа завершён, когда достигнут целевой угол траектории
        if ~ds.ph2_done && st.gamma <= p.theta_target
            ds.ph2_done = true;
        end
    else
        ds.phase = 1;
    end

    % --- Целевой угол атаки для текущей фазы ---
    switch ds.phase
        case 2
            if ds.ph2_done
                alpha_target = 0.0;
            else
                alpha_target = p.alpha_ph2;
            end
        case 3
            alpha_target = p.alpha_apogee;
        case 4
            alpha_target = p.alpha_low;
        otherwise
            alpha_target = 0.0;
    end

    % --- Запуск сглаженной перекладки при смене цели ---
    if abs(alpha_target - ds.alpha_cmd) > 1e-9
        % Текущее фактическое значение становится точкой отсчёта
        s_now = smoothstep5( (tk - ds.t_switch)/p.t_alpha_smooth );
        ds.alpha_from = ds.alpha_from + (ds.alpha_cmd - ds.alpha_from)*s_now;
        ds.alpha_cmd  = alpha_target;
        ds.t_switch   = tk;
    end
end
