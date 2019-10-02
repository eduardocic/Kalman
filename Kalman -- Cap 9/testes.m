clear all; close all; clc;

% F = [0 1 0 0;
%      0 0 0 0;
%      0 0 0 1;
%      0 0 0 0]; 
syms phi t us

% Determinação do parâmetros da matriz 'Q'.
a1 = [0;
      us;
      0;
      us];
a2 = [0 us 0 us];
Q  = a1*a2

% Efetua as multiplicações necessárias.
Q = [0 0 0 0;
     0 1 0 1;
     0 0 0 0;
     0 1 0 1];
 
F = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 0 0]; 
 
PHI  = eye(4) + F*t
PHIt = eye(4) + F'*t   % Cálculo da transposta.


% Cálculo -- Q{k}
PHI * Q * (PHIt)