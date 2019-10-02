function PHI_k = determina_matriz_PHI_k (Ts, f21, f22)

% Primeira linha.
PHI_k(1,1) = 1;
PHI_k(1,2) = Ts;

% Segunda linha.
PHI_k(2,1) = f21*Ts;
PHI_k(2,2) = (1 + f22*Ts);

end