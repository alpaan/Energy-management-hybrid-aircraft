function df_dm = df_m(m, M, N, alpha, beta, alpha_gen)
%DF_PHI_PARALLEL Summary of this function goes here
%   Detailed explanation goes here
df_dm = beta(2)* ( 2 *M* diag(m) + N ); 
end

