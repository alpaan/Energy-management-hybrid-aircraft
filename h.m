function P_c = h(P_em, alpha)
% Function: electric motor loss map between the em output power P_em (MW) 
% to em input power P_c (MW) 
% Input: 
%         - P_em: electric motor output power P_em (MW)
%         - alpha: loss map coefficients 
% Output: 
%         - P_c: electric motor input power P_c (MW) 
out = 0; 
if (alpha(1, 3) == 0)  % linear map
    P_c = alpha(:,2).*P_em + alpha(:,1);
elseif(alpha(1, 3) ~= 0)  % quadratic map
    P_c = alpha(:,3).*P_em.^2 + alpha(:,2).*P_em + alpha(:,1);
end 
end

