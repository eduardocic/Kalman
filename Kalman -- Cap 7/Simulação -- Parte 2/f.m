function xd = f(x, t)

% Constantes utilizadas no cálculo da derivada do estado.
beta = 500;
g    = 32.2;

% Resultado de saída.
xd(1,1) = x(2);
xd(2,1) = 0.0034*g*exp(-x(1)/22000)*(x(2)^2)/(2*beta) - g;

end