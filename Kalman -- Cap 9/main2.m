% =========================================================================
% Feito por: Eduardo H. Santos
%
% Data: 03/10/2019
% =========================================================================
clear; close all; clc;

% Tempo de amostragem do sinal e 
Ts   = 1;
h    = 0.001;

% =========================================================================
%
%       Tratativas realizadas para o objeto perseguido pelo radar
%
% =========================================================================
% Setagem inicial (ângulo e velocidade de disparo).
ang0 = 45;
v0   = 3000;
g    = 32.2;

% Vetor na condição inicial.
X{1} = [0; v0*cos(ang0*pi/180); 0; v0*sin(ang0*pi/180)];


% =========================================================================
%
%                                   Radar
%
% =========================================================================
% Posição do radar.
x_radar = 100000;            
y_radar = 0;

% Característica do ruído do radar relacionado à medida de distância.
sigma_raio  = 100;
media_raio  = 0;

% Característica do ruído do radar relacionado à medida de ângulo.
sigma_theta = 0.01;
media_theta = 0;     



% =========================================================================
%
%                 Informações para o Filtro de Kalman
%
% =========================================================================
% Identidade.
I    = eye(4);

% Matriz 'PHI' (do tipo 'I + F*t')
F    = [0  1  0  0;
        0  0  0  0;
        0  0  0  1;
        0  0  0  0];
% Para o instante de amostragem escolhido:
PHI  = I + (F * Ts);

% Matriz 'R' -- Essa matriz aqui não muda durante a simulação. 
R    = [sigma_theta^2         0;
              0         sigma_raio^2];
          
% Matriz Q{k} -- Aparentemente a equação do livro está errada.
Q = [(Ts^3)/3  (Ts^2)/2      0         0;
     (Ts^2)/2      Ts        0         0;
           0        0    (Ts^3)/3  (Ts^2)/2;
           0        0    (Ts^2)/2     Ts];
PHI_s = 0;
Q     = PHI_s*Q;

% Inicialização do problema (estimativas iniciais para o vetor de estados).
delta    = [ 1000; -100; -1000; 100];
X_hat{1} = X{1} + delta;
        
% Inicialização da matriz de covariância.
P0 = [delta(1)^2     0          0           0;
          0      delta(2)^2     0           0;
          0          0       delta(3)^2     0;
          0          0          0       delta(4)^2];


% =========================================================================
%
%       Roda um Runge-Kuta de 2nd Ordem -- para o sistema sem ruído.
%
% =========================================================================
tempoC(1) = 0;
j  = 0;
i  = 1;
while( X{i}(3,1) >= 0 )
      
    % Base de tempo.
    tempoC(i+1) = i*h;
    
    % Dinâmica do sistema real -- 'X'.
    % --------------------------------
    p1     = f(X{i}, tempoC(i));
    p2     = f(X{i}, tempoC(i) + h);
    X{i+1} = X{i} + (1/2)*h*(p1 + p2);  
    
    % Dinâmica estado estimado -- 'X_hat'.
    % ------------------------------------
    c1         = f(X_hat{i}, tempoC(i));
    c2         = f(X_hat{i}, tempoC(i) + h);
    X_hat{i+1} = X_hat{i} + (1/2)*h*(c1 + c2); 
    
    
    if( rem(tempoC(i+1), Ts) == 0 )
        j         = j + 1;
        tempoD(j) = tempoC(i+1);
        
        % =================================================================
        % Cálculo do 'r' e 'theta' no instante de amostragem -- lido.
        [raio(j), theta(j)] = raio_e_angulo(X{i+1}(1,1), X{i+1}(3,1), x_radar, y_radar);
        
        % Insere o ruído no raio e no ângulo e junta o vetor de medidas.
        theta_noise(j) = theta(j) + (media_theta + sigma_theta*randn);
        raio_noise(j)  =  raio(j) + ( media_raio +  sigma_raio*randn);
        z{j}           = [theta_noise(j); raio_noise(j)];
        
        
        % =================================================================
        % Capturo o estado estimado para o presente instante.
        X_bar{j+1} = X_hat{i+1};
        
        % Apenas quebro as variáveis para o plot em seguida.
        xhat = X_bar{j+1}(1,1);
        yhat = X_bar{j+1}(3,1);

        % Estimação do 'raio' e do 'theta' com base nas medidas até então
        % disponíveis.
        [raio_estimado(j), theta_estimado(j)] = raio_e_angulo(xhat, yhat, x_radar, y_radar);
        % =================================================================
        
        % O meu problema de exige uma matriz 'H', dada por:
        H   = matrizH(xhat, yhat, x_radar, y_radar);
        
        % Resíduo.
        Residuo{j} = z{j} - [theta_estimado(j); raio_estimado(j)];        

        % Equações de Riccati
        if( j == 1 )
            M{j} = PHI * P0 * (PHI') + Q;
            K{j} = M{j} * (H') * inv((H * M{j} * (H') + R));
            P{j} = (I - K{j} * H) * M{j};
        else
            M{j} = PHI * P{j-1} * (PHI') + Q;
            K{j} = M{j} * (H') * inv((H * M{j} * (H') + R));
            P{j} = (I - K{j} * H) * M{j};
        end
         
        % Atualiza o vetor de estados estimados.
        X_hat{i+1} = X_bar{j+1} + K{j}*Residuo{j};       
        
        
    end
    
    % Atualiza o passo de integração.
    i = i + 1;
end

% Em 'x_hat'
for i = 1:(max(size(X_hat)))
    xh(i)  = X_hat{i}(1,1);
    vxh(i) = X_hat{i}(2,1);
    yh(i)  = X_hat{i}(3,1);
    vyh(i) = X_hat{i}(4,1);
end

% Em 'x'
cont = 0;
for i = 1:(max(size(X)))
%     if ( rem(tempoC(i+1), Ts) == 0 )
        cont = cont + 1;
        x(cont)  = X{i}(1,1);
        vx(cont) = X{i}(2,1);
        y(cont)  = X{i}(3,1);
        vy(cont) = X{i}(4,1); 
%     end
end

% De forma equivalente, os valores das covariâncias do sistema.
for k = 1:max(size(P))
      p11(k) = P{k}(1,1);
      p22(k) = P{k}(2,2);
      p33(k) = P{k}(3,3);
      p44(k) = P{k}(4,4);
end

% Pela figura acima, não dá para se ter um pleno entendimento do 
% funcionamento do Filtro de Kalman Estendido. Para melhor entender o 
% resultado obtido, vamos fazer um plot paralelo das variáveis.
% Erro em posição 'x'
figure;
erro_x = x - xh;
subplot(2,2,1);
plot(tempoC, erro_x); 
hold on;
plot(tempoD, sqrt(p11), '--k');
plot(tempoD, -sqrt(p11), '--k');
grid;
xlabel('Tempo (s)');
ylabel('Erro em posição -- x');
ylim([-170 170]);


% Erro em posição 'y'
erro_y = y - yh;
subplot(2,2,2);
plot(tempoC, erro_y);
hold on;
plot(tempoD, sqrt(p33), '--k');
plot(tempoD, -sqrt(p33), '--k');
grid;
xlabel('Tempo (s)');
ylabel('Erro em posição -- y');
ylim([-700 700]);


% Erro em velocidade 'vx'
erro_vx = vx - vxh;
subplot(2,2,3);
plot(tempoC, erro_vx);
hold on;
plot(tempoD, sqrt(p22), '--k');
plot(tempoD, -sqrt(p22), '--k');
grid;
xlabel('Tempo (s)');
ylabel('Erro em velocidade -- vx');
ylim([-20 20]);


% Erro em velocidade 'vy'
erro_vy = vy - vyh;
subplot(2,2,4);
plot(tempoC, erro_vy);
hold on;
plot(tempoD, sqrt(p44), '--k');
plot(tempoD, -sqrt(p44), '--k');
grid;
xlabel('Tempo (s)');
ylabel('Erro em velocidade -- vy');
ylim([-20 20]);

% % Salvar os parâmetros do sistema.
% rodada = num2str(PHI_s);
% name   = 'simulacao_PHIs_';
% name   = strcat(name, rodada);
% name   = strcat(name, '.mat');
% save(name, 'X', 'X_hat', 'tempoC', 'tempoD', 'M', 'P', 'K');
