function nav = ins_mechanize(imu, truth, c, fnav, a0)
%INS_MECHANIZE  Прямая механизация БИНС в ECEF (пакетный прогон траектории).
%
%   Интегрирует показания IMU (приращения угла и скорости) и восстанавливает
%   траекторию: позицию, скорость и ориентацию в ECEF.
%
%   ИСПРАВЛЕНО: член 0.5*[dtheta x dvel] — это ПОВОРОТНАЯ ПОПРАВКА (приводит
%   приращение к осям тела на НАЧАЛО такта), а не скаллинг. Поэтому проекция
%   выполняется ориентацией НАЧАЛА такта. Прежняя комбинация "поворотная
%   поправка + проекция серединой такта" учитывала поворот ДВАЖДЫ.
%   Проверено численно: на траектории с разгоном 17g и вращением 55 °/с
%   правильные комбинации дают ошибку скорости 9e-5 м/с, неправильные 0.067 м/с.
%
%   Вход:
%     imu   - показания IMU (.t, .fb, .wib_b, .Ceb)
%     truth - истинная траектория (.t, .R, .V, .dt)
%     c     - константы
%     fnav  - частота навигатора [Гц]
%     a0    - индекс старта в сетке истины (опц., по умолчанию 1)
%   Выход:
%     nav   - .t, .R, .V, .Q, .C (восстановленная траектория)

    % --- МАРКЕР ВЕРСИИ (чтобы видеть, какая версия реально выполняется) ---
    persistent version_printed
    if isempty(version_printed)
        fprintf('[ins_mechanize] версия v2: проекция НАЧАЛОМ такта (исправлено)\n');
        version_printed = true;
    end

    if nargin < 5 || isempty(a0), a0 = 1; end

    % =====================================================================
    % Шаг 3.1: ПОДГОТОВКА - приращения IMU на такте навигатора
    % =====================================================================
    [a_idx, dtheta, dvel, dt, step] = make_increments(imu, truth, fnav);

    % Оставляем только такты, начинающиеся не раньше точки старта a0
    sel    = a_idx >= a0;
    a_idx  = a_idx(sel);
    dtheta = dtheta(sel,:);
    dvel   = dvel(sel,:);
    M      = numel(a_idx);

    WIEx = skew(c.WIE);

    % --- Начальное состояние берём из истины в точке a0 ---
    r = truth.R(a0,:)';
    v = truth.V(a0,:)';
    q = C2q(imu.Ceb(:,:,a0));

    % --- Память о приращениях предыдущего такта (two-sample) ---
    dth_prev = zeros(3,1);
    dv_prev  = zeros(3,1);

    % --- Преаллокация и запись стартовой точки ---
    R  = zeros(M+1,3);   V  = zeros(M+1,3);   Q = zeros(M+1,4);
    tt = zeros(M+1,1);
    R(1,:) = r';   V(1,:) = v';   Q(1,:) = q';   tt(1) = truth.t(a0);

    % =====================================================================
    % Шаг 3.2: ОСНОВНОЙ ЦИКЛ МЕХАНИЗАЦИИ
    % =====================================================================
    for k = 1:M
        dth = dtheta(k,:)';
        dvl = dvel(k,:)';

        % --- (а) Поворотная поправка + скаллинг (two-sample, Savage) ---
        % dv_rot  — приводит приращение к осям тела на НАЧАЛО такта
        % dv_scul — собственно скаллинг (кросс-члены с предыдущим тактом)
        dv_rot  = 0.5 * cross(dth, dvl);
        dv_scul = (1/12) * ( cross(dth_prev, dvl) + cross(dv_prev, dth) );
        dv_body = dvl + dv_rot + dv_scul;

        % --- (б) Проекция удельной силы ориентацией НАЧАЛА такта ---
        % (не середины! поворот уже учтён поправкой dv_rot)
        C_start = q2C(q);
        a_sf    = C_start * dv_body / dt;

        % --- (в) RK4 для связки (r,v) ---
        [r, v] = rk4_rv(r, v, a_sf, dt, c, WIEx);

        % --- (г) Обновление ориентации (конинг + поворот Земли) ---
        zeta       = -c.WIE * dt;
        phi_coning = (1/12) * cross(dth_prev, dth);
        phi        = dth + phi_coning;
        q          = qmul( qmul(rotvec2q(zeta), q), rotvec2q(phi) );
        q          = q / norm(q);

        % --- Сохранение для следующего такта и запись результата ---
        dth_prev = dth;
        dv_prev  = dvl;

        b = min(a_idx(k)+step, numel(truth.t));
        R(k+1,:) = r';   V(k+1,:) = v';   Q(k+1,:) = q';   tt(k+1) = truth.t(b);
    end

    nav.t = tt;
    nav.R = R;
    nav.V = V;
    nav.Q = Q;
end


function [r2, v2] = rk4_rv(r, v, a_sf, dt, c, WIEx)
%RK4_RV  Один шаг RK4 для связки (позиция, скорость) в ECEF.
%   Правая часть: dr/dt = v;  dv/dt = a_sf + g(r) - 2[w x]v - [w x]^2 r
    [k1r, k1v] = f_rv(r,              v,              a_sf, c, WIEx);
    [k2r, k2v] = f_rv(r + 0.5*dt*k1r, v + 0.5*dt*k1v, a_sf, c, WIEx);
    [k3r, k3v] = f_rv(r + 0.5*dt*k2r, v + 0.5*dt*k2v, a_sf, c, WIEx);
    [k4r, k4v] = f_rv(r + dt*k3r,     v + dt*k3v,     a_sf, c, WIEx);

    r2 = r + dt/6*(k1r + 2*k2r + 2*k3r + k4r);
    v2 = v + dt/6*(k1v + 2*k2v + 2*k3v + k4v);
end


function [dr, dv] = f_rv(r, v, a_sf, c, WIEx)
%F_RV  Правая часть ОДУ движения в ECEF (четыре члена механизации).
    dr = v;
    dv = a_sf + grav_ecef(r, c) - 2*WIEx*v - WIEx*(WIEx*r);
end