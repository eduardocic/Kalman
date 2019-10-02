clear all; close all; clc;

%% ========================================================================
%
% 1. Definição inicial das variáveis; e
% 2. Inserção de ruído na variável a ser amostrada.
%
% =========================================================================

% Tempo total de simulação e variáveis de passo de integração.
h      = 0.001;
Ttotal = 30;
Ts     = 0.1;
t      = linspace(0, Ttotal, (Ttotal)/h + 1);  

% Parâmetros iniciais da nossa simulação (altitude e velocidade iniciais).
X{1} = [200000; -6000];

% Coeficiente de arrasto e constante gravitacional.
beta = 500;
g    = 32.2;

% Por meio da equação não-linear da queda de um corpo submetido à gravidade
% e um arrasto com coeficiente de arrasto constante, a gente integra a
% dinâmica e encontraremos a dinâmica da posição e velocidade.
for k = 1:max(size(t))
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
    
    % Pega os valores da estrutura 'x' (para condensar os resultados).
    x(1) = X{k}(1,1);
    x(2) = X{k}(2,1);
      
    % Cálculo de c1.
    x_d    = zeros(2,1);   
      
    x_d(1) = x(2);
    x_d(2) = 0.0034*g*exp(-x(1)/22000)*(x(2)^2)/(2*beta) - g;
    
    c1     = [x_d(1);
              x_d(2)];
        
    % Cálculo de c2 -- Nosso caso não depende do tempo.
    x_d = zeros(2,1);
      
    x_d(1) = x(2);
    x_d(2) = 0.0034*g*exp(-x(1)/22000)*(x(2)^2)/(2*beta) - g; 
    
    c2     = [x_d(1);
              x_d(2)];
    
    % Dinâmica 'final'.
    X{k+1} = X{k} + (1/2)*h*(c1 + c2);
end

% Conforme descrito, para o referido exemplo o ruído de medida apresenta um
% desvio padrão de 1000 e média zero (variável de saída).
desvio_padrao = 25;
ruido         = 0 + desvio_padrao*randn(1 + (Ttotal/Ts),1);

% A medida 'real' do radar (ou seja, a realidade + ruído).
H  = [1 0];
j  = 0;
for k = 1:max(size(t))
   
   if (mod(t(k), Ts) == 0)
       % Atualiza o contador e o tempo discreto.
       j     = j + 1;
       td(j) = (j-1)*Ts;
       
       % Saída amostrada
       z(j)  = H*X{k} + ruido(j); 
   end
   
   % Altitude.
   x(k) = X{k}(1,1);
end

% Plota o sinal real (em vermelho) e o sinal com ruído (azul).
plot(td, z, 'Linewidth', 2); hold on;
plot(t, x, 'r', 'Linewidth', 2); grid;
xlabel('t(s)');
ylabel('Altitude (ft)');
title('\sigma_{v} = 1000ft');
legend('Sinal medido -- com ruído', 'Sinal real -- sem ruído');


%% ========================================================================
%
%                           FILTRO DE KALMAN
%
% =========================================================================
% Inicialmente, vamos falar que existe um 'match' perfeito entre a dinâmica
% real e o nosso modelo matemático. Para nosso primeiro exemplo, faremos:
u_s  = 10;

% Em consequência, tem-se que:
w     = [0; 
         u_s];
PSI_s = u_s^2;

% Função identidade.
I = eye(2);

% O ruído de medida é o do próprio radar. No instante inicia k = 1, tem-se:
R{1}  = desvio_padrao^2;

% Inicializando a 'Matriz de Covariância' do nosso sistema.
P{1}(1,1) = desvio_padrao^2;
P{1}(2,2) = 20000;
P0        = P{1};

x1_hat(1)  = 200025;
x2_hat(1) = -6150;

% Vetor de estados estimados no instante inicial (chute inicial).
X_hat{1}  = [x1_hat(1);
             x2_hat(1)];

% Contador para os instantes de amostragem.
for k = 1:max(size(z))

    % Atualiza a escala de tempo.
    t_sample(k) = (k - 1)*Ts;
           
    % Qual é a melhor estimativa que temos do sinal no instante atual?
    X_hat_good_estimate{k} = X_hat{k};
      
    % Quebro a variável em variáveis individuais.
    x1_hat(k) = X_hat_good_estimate{k}(1,1);
    x2_hat(k) = X_hat_good_estimate{k}(2,1);
    
    % Determinação de f21 e f22.
    f21 = -(0.0034 * exp(-x1_hat(k)/22000) * g * x2_hat(k)^2)/ (44000 * beta);
    f22 =  (0.0034 * exp(-x1_hat(k)/22000) * g * x2_hat(k)  )/ (beta);

    % Matriz de erros de estado Q{k}.
    % a) Primeira linha.
    q(1,1) = (Ts^3)/3;
    q(1,2) = (1/2)*(Ts^2) + (1/2) * f22 * (Ts^3);

    % b) Segunda linha.
    q(2,1) = (1/2)*(Ts^2) + (1/3) * f22 * (Ts^3);
    q(2,2) = Ts + f22*(Ts^2) + (1/3) * (f22^2) * (Ts^3);

    % Matriz 'Q{k}'.
    Q{k} = PSI_s*q;

    % Matriz 'R{k}'.
    % ==============
    R{k} = R{1};    

    % Matriz de estados atualizada (para utilização nas eq. de Riccati).
    % ==================================================================
    % a_ Primeira linha.
    phi_k(1,1) = 1;
    phi_k(1,2) = Ts;

    % b) Segunda linha.
    phi_k(2,1) = f21 * Ts;
    phi_k(2,2) = (1 + f22 * Ts);

    % Aproximação de Taylor das equaçõs de estados.
    PHI{k} = phi_k;

    % Determina-se a matriz de Kalman.
    % ================================
    if (k == 1)
        M{k} = PHI{k} * P0 * (PHI{k}') + Q{k}; 
        K{k} = M{k} * (H') * (H * M{k} * (H') + R{k})^(-1);
        P{k} = (I - K{k} * H) * M{k};
    else
        M{k} = PHI{k} * P{k-1} * (PHI{k}') + Q{k}; 
        K{k} = M{k} * (H') * (H * M{k} * (H') + R{k})^(-1);
        P{k} = (I - K{k} * H) * M{k};
    end

    % A gente pode melhorar essa estimativa, por meio do Filtro de Kalman
    % Estentido.
    X_hat_best_estimate{k} = X_hat_good_estimate{k} + K{k}*(z(k) - H*X_hat_good_estimate{k});
    
    % Atualiza o estado do sistema.
    X_hat{k}               = X_hat_best_estimate{k};
    
    % Realiza as integrações como se nada estivesse ocorrendo.
    c1 = f(X_hat{k}, t_sample(k));
    c2 = f(X_hat{k}, t_sample(k) + Ts);   
    X_hat{k+1} = X_hat{k} + (1/2)*Ts*(c1 + c2);
    
end

% Cálculo da medida de erro do sistema.
% =====================================
i = 0;
for k = 1:(max(size(X))-1)
    if(mod(t(k), Ts) == 0)    
        i = i + 1;
        
        erro{i} = X{k} - X_hat{i};
        
        % Erro em distância e velocidade.
        dist(i) = erro{i}(1,1);
        vel(i)  = erro{i}(2,1);
        
        % Matrizes p11 e p22.
        p11(i) = sqrt(P{i}(1,1));
        p22(i) = sqrt(P{i}(2,2));
    end
end

% Ajuste e plot do resultado.
% ===========================
% a) Erro em posição.
figure;
plot(t_sample, dist); 
hold on;
plot(t_sample, p11, '--k');
plot(t_sample, -p11, '--k');
grid;
xlabel('t(s)');
ylabel('Erro em altitude (ft)');
title('\Phi_s = 0 e \sigma_{v} = 1000ft');
legend('Erro em posição', 'P_{11}^{1/2}');

% b) Erro em velocidade.
figure;
plot(t_sample, vel);
hold on;
plot(t_sample, p22, '--k');
plot(t_sample, -p22, '--k');
grid;
xlabel('t(s)');
ylabel('Erro em velocidade (ft/s)');
title('\Phi_s = 0 e \sigma_{v} = 1000ft');
legend('Erro em velocidade', 'P_{22}^{1/2}');