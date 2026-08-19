function C = C_e_n(lat, lon)
%C_E_N  Матрица перехода из локальной ENU в ECEF.
%   Столбцы матрицы - орты локального трёхгранника [Восток, Север, Вверх],
%   выраженные в осях ECEF. Применение: r_ecef = C * r_enu.
%
%   Вход:  lat, lon - широта, долгота точки [рад]
%   Выход: C (3x3) - матрица ENU -> ECEF
    sL = sin(lat);  cL = cos(lat);
    sO = sin(lon);  cO = cos(lon);

    East  = [-sO;     cO;     0 ];   % орт "Восток"  в ECEF
    North = [-sL*cO; -sL*sO;  cL];   % орт "Север"   в ECEF
    Up    = [ cL*cO;  cL*sO;  sL];   % орт "Вверх"   в ECEF

    C = [East, North, Up];
end
