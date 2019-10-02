clear all; close all; clc;

% h ........ Passo de simulação;
% Ts ....... Tempo de amostragem (sampling); e
% Ttotal ... Tempo total de simulação.
h      = 0.001;
Ts     = 0.1;
Ttotal = 30;

% Parâmetros iniciais da nossa simulação (altitude e velocidade iniciais).
X{1} = [200000; -6000];

% Coeficiente de arrasto e constante gravitacional.
beta = 500;
g    = 32.2;

% Conforme descrito, para o referido exemplo o ruído de medida apresenta um
% desvio padrão de 1000 e média zero (variável de saída).
desvio_padrao = 1000;
% desvio_padrao = 25;
ruido         = 0 + desvio_padrao*randn(1 + (Ttotal/h),1);

% Por meio da equação não-linear da queda de um corpo submetido à gravidade
% e um arrasto com coeficiente de arrasto constante, a gente integra a
% dinâmica e encontraremos a dinâmica da posição e velocidade.
for k = 1:max(size(ruido))
    % ==============================================================
    % As equações de Runge-Kuta de 2º ordem são simplificarmente, 
    % expressas por meio de:
    %
    %     c1 = f(x, t);
    %     c2 = f(x, t + h);
    %    
    %     x(k+1) = x(k) + (1/2)*h*(c1 + c2);
    % 
    % O que a gente tentará fazer aqui é reproduzir isso no inline.
    % ==============================================================
    % Escala de tempo.
    t(k) = (k-1)*h;
    
    % Coeficientes.
    c1 = f(X{k}, t(k));
    c2 = f(X{k}, t(k) + h);
    
    % Atualiza o estado para 'k+1'.
    X{k+1} = X{k} + (1/2)*h*(c1 + c2);
end

% A medida 'real' do radar (ou seja, a realidade + ruído).
C = [1 0];
for k = 1:max(size(t)) 
   % Saída.
   y(k) = C*X{k} + ruido(k); 
   
   % Altitude.
   x(k) = X{k}(1,1);
end

% % Plota o sinal real (em vermelho) e o sinal com ruído (azul).
% plot(t, y, 'Linewidth', 2); hold on;
% plot(t, x, 'r', 'Linewidth', 2); grid;
% xlabel('t(s)');
% ylabel('Altitude (ft)');
% title('\sigma_{v} = 1000ft');
% legend('Sinal medido -- com ruído', 'Sinal real -- sem ruído');

% =========================================================================
%
%                           FILTRO DE KALMAN
%
% =========================================================================
% Inicialmente, vamos falar que existe um 'match' perfeito entre a dinâmica
% real e o nosso modelo matemático. Para nosso primeiro exemplo, faremos:
u_s  = 0;

% Em consequência, tem-se que:
w     = [0; 
         u_s];
PSI_s = u_s^2;

% O ruído de medida é o do próprio radar. No instante inicia k = 1, tem-se:
R{1}  = desvio_padrao^2;

% Inicializando a 'Matriz de Covariância' do nosso sistema.
P{1}(1,1) = desvio_padrao^2;
P{1}(2,2) = 20000;
P0        = P{1};

x1_hat(1) = 200025;
x2_hat(1) = -6150;


% Vetor de estados no instante inicial.
X_hat{1}  = [x1_hat(1);
             x2_hat(1)];

cont = 0;         
for k = 1:max(size(ruido))
    
    % Atualiza o instante de tempo.
    t(k) = (k-1)*h;
    
    % ==============================================================
    % As equações de Runge-Kuta de 2º ordem são simplificarmente, 
    % expressas por meio de:
    %
    %     c1 = f(x, t);
    %     c2 = f(x, t + h);
    %    
    %     x(k+1) = x(k) + (1/2)*(c1 + c2);
    % 
    % O que a gente tentará fazer aqui é reproduzir isso no inline.
    % ==============================================================
    c1 = f(X_hat{k}, t(k));
    c2 = f(X_hat{k}, t(k) + h);
    X_hat{k+1} = X_hat{k} + (1/2)*h*(c1 + c2);
    
    % Nos instantes de amostragem (Ts)
    if (mod(t(k), Ts) == 0)    
        % Atualiza o contador.
        cont = cont + 1;
        
        % Se nada ocorresse, no instante 'k+1' essa seria a nossa melhor
        % estimativa do sinal baseado apenas nas melhores estimativas
        % disponíveis até o instante 'k'. Sendo assim, pego essas
        % informações.
        X_good{cont} = X_hat{k+1};

        % Além disso, eu preciso atualizar as informações de 'f21' e 'f22'
        % com as melhores estimativas disponíveis até o instante 'k'.
        % Sendo assim, tem-se que:
        x1_hat(cont) = X_hat{k}(1,1);
        x2_hat(cont) = X_hat{k}(2,1);
        
        % E portanto:
        f21 = -(0.0034 * exp(-x1_hat(cont)/22000) * g * x2_hat(cont)^2)/ (44000 * beta);
        f22 =  (0.0034 * exp(-x1_hat(cont)/22000) * g * x2_hat(cont)  )/ (beta);
        
        % Matriz de erros de estado Q{k}.
        % a) Primeira linha.
        q(1,1) = (Ts^3)/3;
        q(1,2) = (1/2)*(Ts^2) + (1/2) * f22 * (Ts^3);

        % b) Segunda linha.
        q(2,1) = (1/2)*(Ts^2) + (1/3) * f22 * (Ts^3);
        q(2,2) = Ts + f22*(Ts^2) + (1/3) * (f22^2) * (Ts^3);

        % Matriz 'Q{k}'.
        Q{cont} = PSI_s*q;

        % Matriz 'R{k}'.
        R{cont} = R{1};    

        % Matriz de estados atualizada (para utilização nas eq. de Riccati).
        % a_ Primeira linha.
        phi_k(1,1) = 1;
        phi_k(1,2) = Ts;

        % b) Segunda linha.
        phi_k(2,1) = f21 * Ts;
        phi_k(2,2) = (1 + f22 * Ts);

        % Aproximação de Taylor das equaçõs de estados.
        PHI{cont} = phi_k;

        % Determina-se a matriz de Kalman.
        if (cont == 1)
            M{cont} = PHI{cont} * P0 * (PHI{cont}') + Q{cont}; 
            K{cont} = M{cont} * (C') * (C * M{cont} * (C') + R{cont})^(-1);
            P{cont} = (eye(2) - K{cont} * C) * M{cont};
        else
            M{cont} = PHI{cont} * P{cont-1} * (PHI{cont}') + Q{cont}; 
            K{cont} = M{cont} * (C') * (C * M{cont} * (C') + R{cont})^(-1);
            P{cont} = (eye(2) - K{cont} * C) * M{cont};
        end
        
        % A gente pode melhorar essa estimativa, por meio do Filtro de Kalman
        % Estentido.
        X_best{cont} = X_good{cont} + K{cont}*(y(k) - C*X_good{cont});

        % Essa melhor estimativas do estado após a passagem do Filtro de
        % Kalman
        % =====================================================================
        X_hat{k+1}  = X_best{cont};
    end
end

% Cálculo da medida de erro do sistema.
% =====================================
for k = 1:(max(size(X))-1)
    erro{k} = X{k} - X_hat{k};
    dist(k) = erro{k}(1,1);
    vel(k)  = erro{k}(2,1);
end

% Ajuste e plot do resultado.
% ===========================
figure;
subplot(2,1,1);
plot(t, dist); grid;
xlabel('t(s)');
ylabel('Erro em altitude (ft)');
title('\Phi_s = 0 e \sigma_{v} = 1000ft');

subplot(2,1,2);
plot(t, vel); grid;
xlabel('t(s)');
ylabel('Erro em velocidade (ft/s)');
title('\Phi_s = 0 e \sigma_{v} = 1000ft');