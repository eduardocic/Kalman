% =========================================================================
% Feito por: Eduardo H. Santos
%
% Data: 25/09/2019
% =========================================================================
clear; close all; clc;




% Parâmetros para a simulação.
% ============================
Tsim = 130;
Ts   = 1;
h    = 0.001;




% Setagem inicial (ângulo e velocidade de disparo).
% =================================================
ang0 = 45;
v0   = 3000;
g    = 32.2;

% Condições iniciais.
x0  = 0;
vx0 = v0*cos(ang0*pi/180);
y0  = 0;
vy0 = v0*sin(ang0*pi/180);

% Vetor na condição inicial.
X{1} = [x0; vx0; y0; vy0];





% Roda um Runge-Kuta de 2nd Ordem -- para o sistema sem ruído.
% ============================================================
for i = 1:((Tsim/h) + 1)
    t(i) = (i-1)*h;
    p1   = f(X{i}, t(i) + h);
    p2   = f(X{i}, t(i));
    
    X{i+1} = X{i} + (1/2)*h*(p1 + p2);
    xt(i)  = X{i}(1,1);
    yt(i)  = X{i}(3,1);
end

% % Plota o resultado -- apenas para conferência.
% % ---------------------------------------------
% plot(xt, yt);
% xlabel('Pos x');
% ylabel('Pos y');
% grid;





% Tratativas com o radar por perto.
% =================================
% Posição do radar.
xr = 100000;            
yr = 0;

for i = 1:max(size(xt))
    % Cálculo do 'r'.
    r(i)     = sqrt((yt(i) - yr)^2 + (xt(i) - xr)^2);
    
    % Cálculo em cima do 'theta', uma vez que a função 'atan2' não
    % normaliza o angulo no intervalor [0, 2*pi] e sim em [-2*pi, 2*pi].
    theta(i) = atan2((yt(i) - yr), (xt(i) - xr));
    if (theta(i) < 0)
        theta(i) = 2*pi + theta(i);
    end
    
    px(i)    = r(i)*cos(theta(i)) + xr;
    py(i)    = r(i)*sin(theta(i)) + yr;
end

% % Plota o resultado -- apenas para conferência.
% % ---------------------------------------------
% figure;
% plot(px, py);
% xlabel('Pos. x');
% ylabel('Pos. y');
% grid;



% Inserção do ruído nas medidas de ângulo e de distância
% ======================================================
sigma_r     = 100;
sigma_theta = 0.01;
media_r     = 0;
media_theta = 0;

% Quantidade de amostras -- n.
QuantSamples = ((Tsim/Ts) + 1);       

ruido_r     = (    media_r + sigma_r*rand(QuantSamples, 1))';
ruido_theta = (media_theta + sigma_theta*rand(QuantSamples, 1))';

% Nos instantes de amostragem eu faço a inserção do ruído na medida real.
j = 0;      % inicializo o contador.
for i = 1:max(size(xt))
    if (mod(t(i), Ts) == 0)
        % Atualizo o contador.
        j = j + 1;
        
        % Insiro um ruído.
        r_ruidoso(j)     =     r(i) + ruido_r(j);
        theta_ruidoso(j) = theta(i) + ruido_theta(j);
        
        % Reconstruo o sinal em cima das medidas ruidosas.
        px_ruidoso(j) = r_ruidoso(j)*cos(theta_ruidoso(j)) + xr;
        py_ruidoso(j) = r_ruidoso(j)*sin(theta_ruidoso(j)) + yr;    
    end
end

% % Plota resultado -- apenas para conferência.
% % -------------------------------------------
% figure;
% plot(px, py, px_ruidoso, py_ruidoso, 'r');
% xlabel('Dist. x');
% ylabel('Dist. y');
% grid;
% legend('Sinal SEM ruído (original)', 'Sinal COM ruído (amostrado pelo radar)');




% Determinação de algumas matrizes do sistema para Kalman Filter.
% ===============================================================
% Matriz 'PHI' (do tipo 'I + F*t')
F    = [0  1  0  0;
        0  0  0  0;
        0  0  0  1;
        0  0  0  0];
% Para o instante de amostragem escolhido, tem-se que:
I    = eye(4);
PHI  = I + (F * Ts);

% Matriz 'R' -- Essa matriz aqui não muda durante a simulação. 
R    = [sigma_theta^2     0;
              0         sigma_r^2];
          
% Matriz Q{k} -- Aparentemente a equação do livro está errada.
% Esta que segue de algumas deduções (considere ela correta).
% Q = [(Ts^3)/3  (Ts^2)/2  (Ts^3)/3  (Ts^2)/2;
%      (Ts^2)/2      Ts    (Ts^2)/2     Ts;
%      (Ts^3)/3  (Ts^2)/2  (Ts^3)/3  (Ts^2)/2;
%      (Ts^2)/2      Ts    (Ts^2)/2     Ts];

Q = [(Ts^3)/3  (Ts^2)/2      0         0;
     (Ts^2)/2      Ts        0         0;
           0        0    (Ts^3)/3  (Ts^2)/2;
           0        0    (Ts^2)/2     Ts];
PHI_s = 0;
Q     = PHI_s*Q;

% Inicialização do problema (estimativas iniciais para o vetor de estados).
delta    = [ 1000; -100; -1000; 100];
x_hat{1} = X{1} + delta;
        
% Inicialização da matriz de covariância.
P0 = [delta(1)^2     0          0           0;
          0      delta(2)^2     0           0;
          0          0       delta(3)^2     0;
          0          0          0       delta(4)^2];
 
cont = 0;     % Contador
for k = 1:((Tsim/h) + 1)
% for k = 1:1
   
    % Base de tempo.
    tempo(k) = (k-1)*h;
    
    % Nos instantes de amostragem, temos de fazer o nosso ajuste baseado:
    % 
    % 1. Medidas obtidas pelo radar (com ruído); e
    % 2. Pelo modelo matemático.
    % if ((mod(tempo(k), Ts) == 0) & (tempo(k) > 0))
    if ((mod(tempo(k), Ts) == 0) & (tempo(k) > 0))
        
        % Atualizo o contador.
        cont = cont + 1;
        tt(cont) = tempo(k);
        
        % Qual a melhor estimativa do vetor de estados até o presente
        % momento? Isso mesmo: 'x_hat'.
        x_good{cont} = x_hat{k};
        
        % Apenas quebro as variáveis para o plot em seguida.
        x_de_t_hat(cont) = x_good{cont}(1,1);
        y_de_t_hat(cont) = x_good{cont}(3,1);
        
        % Eu consigo reconstruir o 'theta' e o 'r' por meio da melhor
        % estimativa e da física do problema, por meio de:
        r_good(cont)     = sqrt((y_de_t_hat(cont) - yr)^2 + (x_de_t_hat(cont) - xr)^2);
        
        % Cálculo em cima do 'theta', uma vez que a função 'atan2' não
        % normaliza o angulo no intervalor [0, 2*pi] e sim em [-2*pi, 2*pi]
        theta_good(cont) = atan2((y_de_t_hat(cont) - yr), (x_de_t_hat(cont) - xr));
        if (theta_good(cont) < 0)
           theta_good(cont) = 2*pi + theta_good(cont);
        end
       
        % O meu problema de exige uma matriz 'H', dada por:
        H   = matrizH(x_de_t_hat(cont), y_de_t_hat(cont), xr, yr);
      
        % Equações de Riccati
        if( cont == 1 )
            M{cont} = PHI * P0 * (PHI') + Q;
            K{cont} = M{cont} * (H') * inv((H * M{cont} * (H') + R));
            P{cont} = (I - K{cont} * H) * M{cont};
        else
            M{cont} = PHI * P{cont-1} * (PHI') + Q;
            K{cont} = M{cont} * (H') * inv((H * M{cont} * (H') + R));
            P{cont} = (I - K{cont} * H) * M{cont};
        end
        
        % Vamos melhorar a estimativa.
        diferenca{cont}    = [(theta_ruidoso(cont) - theta_good(cont));
                              (r_ruidoso(cont)     - r_good(cont))];

        % Esta parte é muito importante, pois a gente utilizará:
        %
        % 1. A melhor estimativa disponível até então (x_good);
        % 2. O cálculo dado o modelo matemático disponível
        x_hat{k}     = x_good{cont} + K{cont} * diferenca{cont};
        x_melhor_estimado{cont} = x_hat{k};
    end
    
    % Se não for o instante de amostragem, continuem a utilizar a
    % integração, utilizando as estimativas do sinal.
    p1   = f(x_hat{k}, t(k));
    p2   = f(x_hat{k}, t(k) + h);
    x_hat{k+1} = x_hat{k} + (1/2)*h*(p1 + p2);
end
    
    


% Pego as variáveis do sistema (já corrigidos).
for i = 1:max(size(x_melhor_estimado))
      xh(i)  = x_melhor_estimado{i}(1,1);
      vxh(i) = x_melhor_estimado{i}(2,1);
      yh(i)  = x_melhor_estimado{i}(3,1);      
      vyh(i) = x_melhor_estimado{i}(4,1);
end


% Pego do estado perfeito as variáveis de contato.
contador = 0;
for i = 1:max(size(xt))
    % Nos instantes de amostragem, temos de fazer o nosso ajuste baseado:
    % 
    % 1. Medidas obtidas pelo radar (com ruído); e
    % 2. Pelo modelo matemático.
    if (mod(tempo(i), Ts) == 0)
        contador = contador + 1;
        xt_a(contador)  = X{i}(1,1);
        vx_a(contador)  = X{i}(2,1);
        yt_a(contador)  = X{i}(3,1);
        vy_a(contador)  = X{i}(4,1);
    end
end



% De forma equivalente, os valores das covariâncias do sistema.
for k = 1:max(size(P))
      p11(k) = P{k}(1,1);
      p22(k) = P{k}(2,2);
      p33(k) = P{k}(3,3);
      p44(k) = P{k}(4,4);
end

% % Faço o plot das variáveis em termos de valores reais e também em termos dos valores estimados.
% figure;
% plot(px, py, xh, yh, 'r');
% xlabel('Dist. x');
% ylabel('Dist. y');
% grid;
% legend('Sinal SEM ruído (original)', 'Sinal estimado (Filtro de Kalman Estendido)');
% title('Construção por coordenadas Cartesianas');

% Pela figura acima, não dá para se ter um pleno entendimento do funcionamento do Filtro de Kalman Estendido. Para melhor entender o resultado obtido, vamos fazer um plot paralelo das variáveis.
% Erro em posição 'x'

xt_a = xt_a(1:(end-1));
yt_a = yt_a(1:(end-1));
vx_a = vx_a(1:(end-1));
vy_a = vy_a(1:(end-1));




figure;
erro_x = xt_a - xh;
subplot(2,2,1);
plot(tt, erro_x);
grid;
xlabel('Tempo (s)');
ylabel('Erro em posição -- x');

% Erro em posição 'y'
erro_y = yt_a - yh;
subplot(2,2,2);
plot(tt, erro_y);
grid;
xlabel('Tempo (s)');
ylabel('Erro em posição -- y');

% Erro em velocidade 'vx'
erro_vx = vx_a - vxh;
subplot(2,2,3);
plot(tt, erro_vx);
grid;
xlabel('Tempo (s)');
ylabel('Erro em velocidade -- vx');

% Erro em velocidade 'vy'
erro_vy = vy_a - vyh;
subplot(2,2,4);
plot(tt, erro_vy);
grid;
xlabel('Tempo (s)');
ylabel('Erro em velocidade -- vy');

