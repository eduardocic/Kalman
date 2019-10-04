function [raio, angulo] = raio_e_angulo(x, y, x_radar, y_radar)

raio   = sqrt((y - y_radar)^2 + (x - x_radar)^2);

angulo = atan2((y - y_radar), (x - x_radar));
if (angulo < 0)
    angulo = 2*pi + angulo;
end

end