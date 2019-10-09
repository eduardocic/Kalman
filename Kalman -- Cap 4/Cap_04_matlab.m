%% Kalman Filter -- Cap 4
% A ideia desse reporte é apresentar o caso do Filtro de Kalman. Estou seguindo 
% a ideia apregoada no livro do Zarchan (_Fundamentals of Kalman Filtering: A 
% Pratical Approach_). Nisso eu irei direto ao problema no tempo contínuo e depois 
% simplifico para o caso em que a gente discretiza o sistema.
% 
% Seja, então, a descrição da física do problema espressa na Espaço de Estados 
% a tempo contínuo:
% 
% $$\dot{x} = Ax + Bu + w$$
% 
% em que $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mi 
% mathvariant="italic">w</mi></mrow></math>$ é um ruído branco do processo. Existe 
% uma matriz associada ao ruído $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi mathvariant="italic">w</mi></mrow></math>$, tal que:
% 
% $$Q = E[ww^T]$$
% 
% Não existe um significado físico para essa matriz, mas a gente verá no 
% futuro que ela está relacionada ao desconhecimento que temos a respeito do processo 
% em si (afinal, não há como representar com 100% de certeza uma dinâmica). 
% 
% O Filtro de Kalman necessita que tenhamos as medidas físicas relacionadas 
% com a seguinte equação:
% 
% $$y = Cx + v$$
% 
% em que $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mi 
% mathvariant="italic">y</mi></mrow></math>$ diz respeito ao vetor de medidas, 
% $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mi 
% mathvariant="italic">C</mi></mrow></math>$é a matriz de medidas e $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi mathvariant="italic">v</mi></mrow></math>$é o ruído 
% branco de medida. Existe uma matriz associada ao ruído $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi mathvariant="italic">v</mi></mrow></math>$, tal que:
% 
% $$R = E[vv^T]$$
% 
% Para utilizarmos o Filtro de Kalman, precisamos discretizar esse sistema 
% expresso a tempo contínuo. Para o sistema, tem-se que:
% 
% $$\Phi(t) = \mathcal{L}^{-1}[(sI-A)^{-1}]$$
% 
% Discretizando a matriz $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi>Φ</mi><mtext> </mtext><mrow><mo>(</mo><mrow><mi mathvariant="italic">t</mi></mrow><mo>)</mo></mrow></mrow></math>$ 
% sistema, tem-se que:
% 
% $$\Phi_{k} = \Phi(T_{s})$$
% 
% em que $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi 
% mathvariant="italic">T</mi></mrow><mrow><mi mathvariant="italic">s</mi></mrow></msub></mrow></math>$ 
% é o período de amostragem. Além do mais, na discretização da nossa amotra, tem-se 
% que:
% 
% $$y_{k} = Cx_k + v_k$$
% 
% e
% 
% $$R_k = E[v_kv_k^T]$$
% 
% A resultante do Filtro de Kalman é dado por:
% 
% $$\hat{x}_{k} = \Phi_k \hat{x}_{k-1} + B_{k}u_{k-1} + K_{k}[y_{k} - H(\Phi_{k}\hat{x}_{k-1}+B_{k}u_{k-1})]$$
% 
% em que $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi 
% mathvariant="italic">u</mi></mrow><mrow><mi mathvariant="italic">k</mi><mo>−</mo><mn>1</mn></mrow></msub></mrow></math>$ 
% é assumido constante entre as amostras, $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><msub><mrow><mi mathvariant="italic">K</mi></mrow><mrow><mi 
% mathvariant="italic">k</mi></mrow></msub></mrow></math>$ é o ganho de Kalman 
% no instante $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mi 
% mathvariant="italic">k</mi></mrow></math>$ e a matriz $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><msub><mrow><mi mathvariant="italic">B</mi></mrow><mrow><mi 
% mathvariant="italic">k</mi></mrow></msub></mrow></math>$ é dado por:
% 
% $$B_{k} = \int_{0}^{T_s}\Phi({\tau})\textbf{B}d\tau$$
% 
% Os ganhos de Kalman são computados por meio da matriz de Riccati. As equações 
% são recursivas:
% 
% $$M_k = \Phi_{k}P_{k-1}\Phi_{k}^T + Q_k$$
% 
% $$K_k = M_{k}C^T(CM_kC^T+R_k)^{-1}$$
% 
% $$P_k = (I-K_kC)M_k$$
% 
% em que:
% 
% $$Q_k = \int_{0}^{T_s}\Phi(\tau)\textbf{Q}\Phi^{T}(\tau) d\tau$$
% 
% Para iniciar as equações de Riccati a gente precisa inicializar a matriz 
% com algum $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi 
% mathvariant="italic">P</mi></mrow><mrow><mn>0</mn></mrow></msub></mrow></math>$.
% 
% Olhando de maneira estrita, o que a gente tem é que:
% 
% $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mo>•</mo></mrow></math>$a 
% matriz* *$<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi 
% mathvariant="italic">M</mi></mrow><mrow><mi mathvariant="italic">k</mi></mrow></msub></mrow></math>$ 
% diz respeito à matriz de Covariância representando os erros das estimativas 
% dos estados antes da atualização;
% 
% $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mo>•</mo></mrow></math>$a 
% matriz $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi 
% mathvariant="italic">P</mi></mrow><mrow><mi mathvariant="italic">k</mi></mrow></msub></mrow></math>$ 
% diz respeito à matriz de Covariância representando os erros das estimativas 
% dos estados depois da atualização.
% 
% O que a equação do Filtro de Kalman faz é simplesmente tomar uma estimativa 
% e corrigir essa estimativa com base nos dados reais amostrados. Sendo assim, 
% tomemos um exemplo.
%% Exemplo:
% Imagine o caso em que desejamos estimar a posição e velocidade de um objeto 
% caindo. Ele inicialmente apresenta altitude de 400.000 $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi mathvariant="normal">ft</mi></mrow></math>$ e velocidade 
% inicial de 6.000 $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mfrac><mrow><mi 
% mathvariant="normal">ft</mi></mrow><mrow><mi mathvariant="italic">s</mi></mrow></mfrac></mrow></math>$ 
% em direção ao um radar posicionado no solo (ou seja, caindo). Neste exemplo 
% não se considera o efeito do arrasto e a aceleração da gravidade é g = 32,2 
% $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mfrac><mrow><mi 
% mathvariant="normal">ft</mi></mrow><mrow><msup><mrow><mi mathvariant="italic">s</mi></mrow><mrow><mn>2</mn></mrow></msup></mrow></mfrac></mrow></math>$ 
% . Vamos assumir que as medidas obtidas pelo radar são ruidosas, apresentando 
% as mesmas um desvio padrão de 1.000 ft com média zero. Neste exemplo considera-se 
% que o radar faz amostragem na taxa de 10Hz (ou seja, 10 amostras em 1 s) durante 
% 30 s. A ideia é construir um Filtro de Kalman para estimar a posição e velocidade 
% do objeto a medida que ele cai e são coletadas informações pelo radar.
% 
% 
% 
% Pela física de Newton, sabe-se que a altura de um objeto que cai com as 
% características apresentadas é dada por:
% 
% $$x = 400000- 6000 - gt^2/2$$
% 
% a velocidade é:
% 
% $$\dot{x} = 6000 - gt$$
% 
% Cadastrando os itens descritos, tem-se:

x0 = 400000;        % altura inicial
v0 = 6000;          % velocidade inicial
g  = 32.2;          % gravidade.
%% 
% Do jeito que foi abordado no início do processo, tem-se:

Ttotal = 30;   % Tempo total de amostragem.
fs     = 10;   % Frequência de amostragem.
t      = linspace(0, Ttotal, (fs*Ttotal + 1));        % Instantes de amostragem.
%% 
% O objeto em queda terá, dada as condições do problema proposto, o seguinte 
% comportamento:

% Posição.
for i = 1:(fs*Ttotal + 1)
    x(i) = x0 - v0*t(i) - 0.5*g*(t(i)^2);        
end
% Velocidade.
v = -v0 - g*t;   

% Plota os resultados de posição e velocidade.
subplot(1, 2, 1);
plot(t, x); grid; 
xlabel('t (s)');
ylabel('Posição (pés)');
axis([0 Ttotal 0 400000]);

subplot(1, 2, 2);
plot(t, v); grid;
xlabel('t (s)');
ylabel('Velocidade (pés/s)');
axis([0 Ttotal -10000 0]);
%% 
% O radar (dispositivo que fará as medições do sistema) apresenta um desvio 
% padrão de 1000 ft com média zero. Dessa forma, as medidas do radar estarão corrompidas 
% com erro (o que será denotado por ruido):

desvio_padrao = 1000;                                 % em pés.
ruido = 0 + desvio_padrao*randn(fs*Ttotal + 1, 1);    % de medida pelo radar.
ruido = ruido';                            
%% 
%  Isso significa dizer que a medida do radar será a medida real corrompida 
% com ruido. Sendo assim, Tem-se que:

% Sinal MEDIDO pelo radar com ruído (posição apenas).
y = x + ruido;
figure;
plot(t, y);
grid;

xlabel('t (s)');
ylabel('Posição (pés)');

%% 
% Vamos agora montar o problema do jeito que preconiza o Filtro de Kalman. 
% 
% O meu sistema no final das contas apresenta uma interessante forma de apresentação, 
% pois tem-se que um objeto caindo sofrerá apenas ação da gravidade:
% 
% $$\ddot{x}(t) = -g$$
% 
% E assim, derivando-se essa aceleração (constante), tem-se que:
% 
% $$<math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow><mrow><mfrac><mrow><msup><mrow><mi 
% mathvariant="normal">d</mi></mrow><mrow><mn>3</mn></mrow></msup></mrow><mrow><mi 
% mathvariant="normal">d</mi><msup><mrow><mi mathvariant="italic">t</mi></mrow><mrow><mn>3</mn></mrow></msup></mrow></mfrac><mrow><mi 
% mathvariant="italic">x</mi><mrow><mo>(</mo><mrow><mi mathvariant="italic">t</mi></mrow><mo>)</mo></mrow></mrow></mrow><mo>=</mo><mn>0</mn></mrow></math>$$
% 
% Então, a gente consegue escrever o sistema no seguinte formato:
% 
% $$<math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow><mrow><mo>[</mo><mtable><mtr><mtd><mrow><mrow><mfrac><mrow><mi 
% mathvariant="normal">d</mi></mrow><mrow><mi mathvariant="normal">d</mi><mrow><mi 
% mathvariant="italic">t</mi></mrow></mrow></mfrac><mrow><mi mathvariant="italic">x</mi></mrow></mrow></mrow></mtd></mtr><mtr><mtd><mrow><mrow><mfrac><mrow><msup><mrow><mi 
% mathvariant="normal">d</mi></mrow><mrow><mn>2</mn></mrow></msup></mrow><mrow><mi 
% mathvariant="normal">d</mi><msup><mrow><mi mathvariant="italic">t</mi></mrow><mrow><mn>2</mn></mrow></msup></mrow></mfrac><mrow><mi 
% mathvariant="italic">x</mi></mrow></mrow></mrow></mtd></mtr><mtr><mtd><mrow><mrow><mfrac><mrow><msup><mrow><mi 
% mathvariant="normal">d</mi></mrow><mrow><mn>3</mn></mrow></msup></mrow><mrow><mi 
% mathvariant="normal">d</mi><msup><mrow><mi mathvariant="italic">t</mi></mrow><mrow><mn>3</mn></mrow></msup></mrow></mfrac><mrow><mi 
% mathvariant="italic">x</mi></mrow></mrow></mrow></mtd></mtr></mtable><mo>]</mo></mrow><mo>=</mo><mrow><mo>[</mo><mtable><mtr><mtd><mrow><mn>0</mn></mrow></mtd><mtd><mrow><mn>1</mn></mrow></mtd><mtd><mrow><mn>0</mn></mrow></mtd></mtr><mtr><mtd><mrow><mn>0</mn></mrow></mtd><mtd><mrow><mn>0</mn></mrow></mtd><mtd><mrow><mn>1</mn></mrow></mtd></mtr><mtr><mtd><mrow><mn>0</mn></mrow></mtd><mtd><mrow><mn>0</mn></mrow></mtd><mtd><mrow><mn>0</mn></mrow></mtd></mtr></mtable><mo>]</mo></mrow><mrow><mo>[</mo><mtable><mtr><mtd><mrow><mi 
% mathvariant="italic">x</mi></mrow></mtd></mtr><mtr><mtd><mrow><mover><mrow><mi 
% mathvariant="italic">x</mi><mtext> </mtext></mrow><mrow><mo>˙</mo></mrow></mover></mrow></mtd></mtr><mtr><mtd><mrow><mover><mrow><mtext> 
% </mtext><mi mathvariant="italic">x</mi></mrow><mrow><mo>¨</mo></mrow></mover></mrow></mtd></mtr></mtable><mo>]</mo></mrow></mrow></math>$$
% 
% Significa dizer então que a matriz B = 0 e a matriz A é a que está logo 
% acima.
% 
% Sendo assim, temos de determinar quem é a matriz $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi>Φ</mi><mrow><mo>(</mo><mrow><mi mathvariant="italic">t</mi></mrow><mo>)</mo></mrow></mrow></math>$ 
% para depois discretizá-la.

% Matriz A do sistema (contínuo).
A = [0 1 0;
     0 0 1;
     0 0 0];
B = [0; 0; 0];

% Calcular a inversa de Laplace.
syms s

I     = eye(3);  % Matriz identidade de ordem 3.
phi_t = ilaplace((s*I - A)^(-1))
%% 
% Quanto a gente discretiza o sistema para um $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><msub><mrow><mi mathvariant="italic">T</mi></mrow><mrow><mi 
% mathvariant="italic">s</mi></mrow></msub></mrow></math>$ conforme especificado, 
% tem-se que:

% Matriz de estados.
Ts    = 1/fs;              % Tempo de amostragem.
phi_k = [1    Ts (Ts^2)/2;
         0     1    Ts;
         0     0    1];
     
% Vocẽ poderia também fazer de forma mais inteligente por meio do comando abaixo:
% phi_k = subs(phi_t, Ts);
%% 
% Vamos agora determinar a matriz Q a tempo discreto. Importante saber quem 
% é a mesma, em tempo contínuo é dada por:
% 
% $$<math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow><mi 
% mathvariant="italic">Q</mi><mo>=</mo><mi mathvariant="italic">E</mi><mrow><mo>[</mo><mrow><msup><mrow><mi 
% mathvariant="normal">ww</mi></mrow><mrow><mi mathvariant="italic">T</mi></mrow></msup></mrow><mo>]</mo></mrow></mrow></math>$$
% 
% Como não há erro de processo envolvido no sistema (ou seja, é considerada 
% essa situação), a nossa matriz $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi mathvariant="italic">Q</mi></mrow></math>$ será nula 
% e quadrada de ordem 3. Quando a gente coloca essa matriz para determinação da 
% matriz $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi 
% mathvariant="italic">Q</mi></mrow><mrow><mi mathvariant="italic">k</mi></mrow></msub></mrow></math>$ 
% dará no final uma matriz nula, quadrada e de ordem 3 (só verificar nas equações). 
% Isso pode ser comprovado substituindo Q nula na equação abaixo:
% 
% $$Q_k = \int_{0}^{T_s}\Phi(\tau)\textbf{Q}\Phi^{T}(\tau) d\tau$$
% 
% Assim, tem-se que:

Q   = zeros(3,3);

% 1) Uma forma de fazer (sendo direto).
Q_k = zeros(3,3);

% 2) Podemos também fazer no 'brute force', Para tal, tem-se que:
% PRODUTO = phi_t * Q * (phi_t');
% Q_k     = int(PRODUTO, [0 Ts]);

%% 
% Olhando o tipo de medida que iremos ter, o radar mede apenas a distância 
% do objeto. Sendo assim, a matrix C do sistema terá apenas uma únia entrada.

C = [1 0 0];
%% 
% O desvio padrão do ruído de medida é dado de tal forma que temos que determinar 
% a matriz $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi 
% mathvariant="italic">R</mi></mrow><mrow><mi mathvariant="italic">k</mi></mrow></msub><mo>.</mo></mrow></math>$ 
% Sendo assim, tem-se que o ruído será equivalente ao erro característico do radar.

R_k = desvio_padrao^2;
%% 
% O nosso sistema está todo resolvido. Precisamos apenas iniciar o processo 
% do Filtro de Kalman e rodar a simulação.

% Inicialização do sistema de Kalman.
P{1}  = 999999999*eye(3);

% Estimativas iniciais dos estados do sistema.
x_hat{1} = [0; 0; 0];
%% 
% Executo a recursão do Filtro de Kalman.

% Recursão de Kalman.
for k = 2:max(size(x))
    
    % Equação de Riccati.
    M{k} = (phi_k) * P{k-1} * (phi_k') + Q_k;
    K{k} = M{k} * (C') * inv(C * M{k} * (C') + R_k);
    P{k} = (I - K{k}*C)*M{k};
    
    % Atualização dos estados.
    x_hat{k} = phi_k*x_hat{k-1} + K{k}*(y(k) - C*phi_k*x_hat{k-1});
end

% Separa os resultados para facilitar o plot.
for i = 1:max(size(x))
   x_end(i)   = x_hat{i}(1,1); 
   xd_end(i)  = x_hat{i}(2,1);
   xdd_end(i) = x_hat{i}(3,1);
end
%% 
% Agora eu ploto todo o sistema.

figure;
plot(t, x_end, t, y, 'r')
legend('Estimado (Kalman)', 'Medido');
xlabel('t (s)');
ylabel('Posição (pés)');
grid;
axis([0 Ttotal 0 400000])

figure;
plot(t, xd_end, t, v, 'r')
legend('Estimado (Kalman)', 'Real');
xlabel('t (s)');
ylabel('Velocidade (pés/s)');
axis([0 Ttotal -10000 0])
grid;

figure;
G2 = -g*linspace(1, 1 ,max(size(x)));
plot(t, xdd_end, t, G2, 'r')
legend('Estimado (Kalman)', 'Real');
xlabel('t (em segundos)');
ylabel('Aceleração (em pés/s)');
ylim([-100 100]);
grid;
%% 
% Perceba que o Filtro de Kalman conseguiu fazer uma boa estimação do meu 
% sistema.