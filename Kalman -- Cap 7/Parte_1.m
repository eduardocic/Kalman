clear all; close all; clc;

% Condição inicial
X{1} = [200000; -6000];
beta = [500, 1000, 2000, 9999999999999];


H      = 0.001;
Ttotal = 30;
Ts     = 0.1;
% Roda o Runge-Kutta de 2º Ordem.
for j = 1:max(size(beta))
    for i = 1:(Ttotal/H)
        t(i) = H*i;
        p1   = fun(beta(j), X{i}, t(i));
        p2   = fun(beta(j), X{i}, t(i) + H);
        X{i+1} = X{i} + 0.5*H*(p1+p2);
    end
    % Separa as variáveis.
    for i = 1:(Ttotal/H)
         x(i)      = X{i}(1,1);
         x_dot(i)  = X{i}(2,1);
    end
    plot(t, x); hold on;
    xlabel('Tempo (s)');
    ylabel('Amplitude');
end
grid;
legend('500', '1000', '2000', 'Infty');
