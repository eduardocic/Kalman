clear all; close all; clc;

% Carrega os arquivos já salvos -- Simulação 1
load('simulacao_PHIs_0.mat');

for k = 1:max(size(P)) 
    p22(k) = P{k}(2,2);
end
sigma_p22_0 = sqrt(p22);


% Carrega os arquivos já salvos -- Simulação 2
load('simulacao_PHIs_10.mat');

for k = 1:max(size(P)) 
    p22(k) = P{k}(2,2);
end
sigma_p22_10 = sqrt(p22);


% Carrega os arquivos já salvos -- Simulação 3
load('simulacao_PHIs_100.mat');

for k = 1:max(size(P)) 
    p22(k) = P{k}(2,2);
end
sigma_p22_100 = sqrt(p22);


% Carrega os arquivos já salvos -- Simulação 4
load('simulacao_PHIs_1000.mat');

for k = 1:max(size(P)) 
    p22(k) = P{k}(2,2);
end
sigma_p22_1000 = sqrt(p22);


plot(tempoD, sigma_p22_0); hold on;
plot(tempoD, sigma_p22_10);
plot(tempoD, sigma_p22_100);
plot(tempoD, sigma_p22_1000);
xlabel('Tempo (s)');
ylabel('Erro na estimativa da velocidade de queda (feet/s)');
legend('\Phi_{s} = 0', '\Phi_{s} = 10', '\Phi_{s} = 100', '\Phi_{s} = 1000');
grid;
