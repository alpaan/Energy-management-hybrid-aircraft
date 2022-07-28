function P_c = inv_f_P_b(P_b, R_U)
% Function: inverse input/output map for the battery electric bus.
% Compute the effective electrical power P_c (MW) delivered to the em for a
% given input battery chemical power P_b (MW).
% Input:  
%         - P_b: battery chemical power (MW) 
%         - R_U: ratio R/U^2*10^6 where R and U are the battery eq.
%         circuit resistance and voltage source respectively 
% Output: 
%         - P_c: effective electrical power (MW)

P_c = P_b - R_U*P_b.^2; 
end

