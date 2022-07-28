function P_b = f_P_b(P_c, R_U)
% Function: input/output relation for the battery electric bus loss map.
% Compute battery chemical power P_b (MW) required to deliver the effective
% electrical power P_c (MW).
% Input:  
%         - P_c: effective electrical power (MW)
%         - R_U: ratio 4*R/U^2*10^6 where R and U are the battery eq.
%         circuit resistance and voltage source respectively 
% Output: 
%         - P_b: battery chemical power (MW)

P_b= 1/(2*R_U)*(1 - sqrt(1 - 4*R_U*P_c)); 
end

