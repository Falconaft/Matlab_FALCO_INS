function r = lla2ecef(lat, lon, h, c)
%LLA2ECEF  Перевод геодезических координат (широта, долгота, высота) в ECEF.
%   Аналог в Navigation Toolbox: lla2ecef (но здесь своя реализация для
%   полного контроля над конвенцией и константами WGS-84).
%
%   Вход:  lat, lon - широта, долгота [рад];  h - высота над эллипсоидом [м]
%   Выход: r - позиция в ECEF [м]
    sL = sin(lat);  cL = cos(lat);
    sO = sin(lon);  cO = cos(lon);

    % Радиус кривизны первого вертикала (prime vertical)
    N = c.A_E / sqrt(1 - c.E2*sL*sL);

    r = [ (N + h)*cL*cO;
          (N + h)*cL*sO;
          (N*(1 - c.E2) + h)*sL ];
end
