clear all; close all; clc;

% Tempo total de simulação e variáveis de passo de integração.
h      = 0.001;
Ttotal = 30;
Ts     = h;
t      = linspace(0, Ttotal, (Ttotal)/h + 1);  

% Parâmetros iniciais da nossa simulação (altitude e velocidade iniciais).
X{1} = [200000; -6000];

% Coeficiente de arrasto e constante gravitacional.
beta = 500;
g    = 32.2;

% Conforme descrito, para o referido exemplo o ruído de medida apresenta um
% desvio padrão de 1000 e média zero (variável de saída).
desvio_padrao = 1000;
ruido         = 0 + desvio_padrao*randn(1 + (Ttotal/h),1);

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
    %     x(k+1) = x(k) + (1/2)*(c1 + c2);
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

% A medida 'real' do radar (ou seja, a realidade + ruído).
C = [1 0];
for k = 1:max(size(t))
   % Saída.
   y(k) = C*X{k} + ruido(k); 
   
   % Altitude.
   x(k) = X{k}(1,1);
end

% Plota o sinal real (em vermelho) e o sinal com ruído (azul).
plot(t, y, 'Linewidth', 2); hold on;
plot(t, x, 'r', 'Linewidth', 2); grid;
xlabel('t(s)');
ylabel('Altitude (ft)');
legend('Sinal medido -- com ruído', 'Sinal real -- sem ruído');


% =========================================================================
%
%                       Filtro de Kalman Extendido
%
% =========================================================================
% 'Ruído' no estado (no caso, apenas na aceleração).

% a) Inicialmente, vamos falar que existe um 'match' perfeito entre a dinâmica
% real e o nosso modelo matemático. 
mu_s  = 0;
w     = [0; mu_s];
PSI_s = mu_s^2;

% b) o ruído de medida é o do próprio radar.
R{1}  = desvio_padrao^2;

% Variáveis auxiliares.
P{1}(1,1) = desvio_padrao^2;
P{1}(2,2) = 20000;

% Melhores estimativas iniciais.
x_hat(1)  = 200025;
xd_hat(1) = -6150;

% Vetor de estados no instante inicial (melhores estimativas).
X_hat{1}  = [x_hat(1);
             xd_hat(1)];

for k = 2:max(size(t))
    
    % Vamos quebrar o problema em partes. Determinares o 'rho' com base nas
    % nossas melhores estimativas. 
    % =====================================================================
    rho(k) = 0.0034*exp(-x_hat(k)/22000);
    
    % Variáveis da matriz de estados (já discretizado). Tais valores serão
    % utilizados para determinação das matrizes a serem utilizadas nos
    % ganhos de Kalman.
    % =====================================================================
    f21(k) = rho(k) * g * (xd_hat(k)^2)/(44000 * beta);
    f22(k) = (rho(k) * (xd_hat(k)) * g)/beta;
    
    % Matriz de erros de estado Q{k}.
    % ===============================
    % a) Primeira linha.
    q(1,1) = (Ts^3)/3;
    q(1,2) = (1/2)*(Ts^2) + (1/2)*f22*(Ts^3);

    % b) Segunda linha.
    q(2,1) = (1/2)*(Ts^2) + (1/3)*f22*(Ts^3);
    q(2,2) = Ts + f22*(Ts^2) + (1/3)*(f22^2)*(Ts^3);

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
    phi_k(2,1) = f21*Ts;
    phi_k(2,2) = (1 + f22*Ts);
    
    % Aproximação de Taylor das equaçõs de estados.
    PHI{i} = phi_k;
   
    
    % Determina-se a matriz de Kalman.
    % ================================
    M{k} = PHI_k{k} * P{k-1} * (PHI{k}') + Q{k}; 
    K{k} = M{k} * (C') * inv(C * M{k} * (C') + R{k});
    P{k} = (eye(2) - K{k} * C) * M{k};
    
    % Dentro dessa perspectiva, a melhor estimativa ATÉ AGORA que a gente
    % pode ter do estado no instante 'k+1' é dada por:
    % =====================================================================
    X_hat_before_kalman{k+1} = X_hat{k} + X_boa_estimativa{k}*Ts;
    
%     % Determina-se a matriz de Kalman.
%     % ===================================
%     [M{i+1}, K{i+1}, P{i+1}] = Kalman(PHI_k{i}, P{i}, H, R_k, Q_k{i});
    
    % Depois dessa determinação da matriz de Kalman, é necessário fazer um
    % ajutes na estimativa que a gente tem no estado estimado. Dessa forma:
    % =====================================================================
    X_hat{k+1} = X_hat_before_kalman{k+1} + K{k+1}*(z(k) - H*X_hat_before_kalman{k+1});
    
    % Desmebra a criança nas variáveis de interesse.
    % ==============================================
    x_hat(k+1)   = X_hat{k+1}(1,1);
    x_d_hat(k+1) = X_hat{k+1}(2,1);
end

% % Pegando os resultados dos erros de estimativa.
% for i = 1:max(size(X))
%     erro{i}     = X{i} - X_hat{i};
%     erro_x(i)   = erro{i}(1,1);
%     erro_d_x(i) = erro{i}(2,1);
% end
% 
% % Erro de posição.
% plot(t, erro_x(1,1:(end-1))); grid;
% figure;
% % Erro de velocidade.
% plot(t, erro_d_x(1,1:(end-1))); grid;
% 
% 
% % Verificação dos ganhos do filtro de Kalman.
% for i = 1:(max(size(K))-1)
%     k1(i) = K{i+1}(1,1);
%     k2(i) = K{i+1}(2,1);
% end
% figure;
% plot(t, k1, t, k2); grid;
% legend('k1', 'k2');
