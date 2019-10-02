function [M_update, K_update, P_update] = Kalman(PHI_k, P_last, H, R_k, Q_k)

M_update = PHI_k*P_last*(PHI_k') + Q_k; 
K_update = M_update*(H')*inv(H*M_update*(H') + R_k);
P_update = (eye(2) - K_update*H)*M_update;

end