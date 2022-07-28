function [f, g, H_inv_vec] = f_eval_P_b(x, z, y, u, mdl)
% Function:  Compute the objective function f(P_b) and (optionally) its 
% gradient and Hessian. 
% Input:
%        - x: optimisation variable (P_b)
%        - z: chi - xi + lambda1
%        - y: drive power
%        - u: zeta + lambda5
%
% Output: 
%         - f: cost function, f(P_b) = 0.5*sigma1*sum((f_phi(m,P_b) + z)^2)
%                                    + 0.5*sigma1*sum((P_b - u)^2)
%         - g: gradient
%         - Hinv_vec: inverse Hessian in vector form

% Penalty parameters (ADMM)
p_chi = mdl.sigma1;
p_zeta = mdl.sigma5; 

% E.m. coefficients
RU = mdl.RU;
b2 = mdl.alpha(:,3); 
b1 = mdl.alpha(:,2); 
b0 = mdl.alpha(:,1);
bem_const = - b1./( 2*b2 ); 

% Fuel coefficients
k1 = mdl.beta(2);
k0 = mdl.beta(1);

% Initialise variables
F_phi = 0;
sqrt_pi1 = 0; 
pi4 = 0; 
f1 = 0; 
f2 = 0; 

% Compute cost function
if (mdl.any == 1)  % if quadratic e.m. loss map 
    pi1 = (-RU.*x.^2 + (x - b0))./b2 + b1.^2./(4*b2.^2);
    sqrt_pi1 = sqrt(pi1);
    pi4 = z + k0 + k1.*(y - bem_const - sqrt_pi1);
    f = sum(p_chi.*(pi4.^2) + p_zeta.*((u - x).^2), 1)/2;
else               % if linear e.m. loss map 
    F_phi = k0 + k1*(y + (RU*x.^2-x + b0)./b1);
    f1 = p_chi/2*(z + F_phi).^2;
    f2 = p_zeta/2*(u-x ).^2; 
    f = sum (f1 + f2); 
end

if nargout > 1  % gradient required
    
    % Compute gradient
    if (mdl.any == 1)  % if quadratic e.m. loss map 
        pi2 = (1 - 2.*RU*x)./b2;
        df_P_b = -k1.*pi2./(2*sqrt_pi1);
        g = (p_zeta.*(x - u)) + p_chi*df_P_b.*pi4;
    else               % if linear e.m. loss map 
        df_P_b = k1*(2*RU*x -1)./b1;
        g = p_chi*df_P_b.*(z + F_phi) -p_zeta*(u - x);
    end

    if nargout > 2 % inverse Hessian vector required 
        
        % Compute Hessian 
        if (mdl.any == 1)  % if quadratic e.m. loss map 
            d2f_P_b = (k1./b2*RU + df_P_b.^2/k1)./sqrt_pi1; 
            H_vec = p_zeta + p_chi*(df_P_b.^2 + d2f_P_b.*pi4); 
        else               % if linear e.m. loss map 
            d2f_P_b = 2*RU*k1./b1; 
            H_vec = p_zeta +  p_chi*d2f_P_b.*(z + F_phi)+ p_chi*(df_P_b.^2);
        end
        
        % Improve conditioning of H before inversion
        tau = max(0, eps^(1/2) - (min(H_vec)));
        H_vec = H_vec + tau;  % make Hessian sufficiently positive definite
        H_inv_vec = 1./H_vec;
    end
end
end