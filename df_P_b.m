function df_dPb = df_P_b(Pb, R, U, alpha, beta, alpha_gen, RU)
%DF_PHI_PARALLEL Summary of this function goes here
%   Detailed explanation goes here 

df_dPb = beta(2)*diag(1./alpha(:,2))* ( 2 * RU * diag(Pb) -eye(length(Pb))); 
end

