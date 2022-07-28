function df_dPb = d2f_P_b(R, U, alpha, beta, alpha_gen, RU)
%DF_PHI_PARALLEL Summary of this function goes here
%   Detailed explanation goes here
df_dPb = 2*RU * beta(2)*diag(1./alpha(:,2)); 
end

