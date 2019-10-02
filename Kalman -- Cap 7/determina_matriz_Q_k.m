function Q_k = determina_matriz_Q_k (PHI_s, Ts, f22)

% Primeira linha.
q(1,1) = (Ts^3)/3;
q(1,2) = 0.5*(Ts^2) + 0.5*f22*(Ts^3);

% Segunda linha.
q(2,1) = 0.5*(Ts^2) + (1/3)*f22*(Ts^3);
q(2,2) = Ts + f22*(Ts^2) + (1/3)*(f22^2)*(Ts^3);

Q_k = PHI_s*q;

end