function [f, g, Hinv_vec] = f_eval_m(x, z, y, u, mdl)
% Function:  Compute the objective function f(m) and (optionally) its 
% gradient and Hessian. 
% Input:
%        - x: optimisation variable (m)
%        - z: chi - xi + lambda1 
%        - y: -g^-1(P_b) where g^-1 is the inverse battery electronics map
%        - u: Psi_xi - m_0_0*Phi + lambda2 
%
% Output: 
%         - f: cost function, f(m) = 0.5*sigma1*sum(f_phi(m, Pb) + z)
%                                  + 0.5*sigma2*sum(m + u)
%                                    
%         - g: gradient
%         - Hinv_vec: inverse Hessian in vector form

% Penalty parameters (ADMM)
p_chi = mdl.sigma1;
p_m = mdl.sigma2;

% Drive power parameters 
eta2 = mdl.m_drv;
eta1 = mdl.n_drv;
eta0 = mdl.q_drv; 

% Fuel map coefficients
k1 = mdl.beta(2);
k0 = mdl.beta(1);


% Compute cost function
pi1 = eta2.*[mdl.m_0_0; x(1:end-1)].^2 + eta1.*[mdl.m_0_0; x(1:end-1)] + eta0 + y;
pi2 = (z + k0 + k1.*pi1);
f = sum(p_chi.*pi2.^2 + p_m.*((x + u).^2), 1)/2;

% Compute radient 
if nargout > 1 % gradient required
    pi3 = 2.*eta2.*[0; x(1:end-1)] + eta1.*[0; ones(size(x(1:end-1)))];
    g = p_m.*(x + u) + p_chi.*pi2.*pi3.*k1;
    
    % Compute inverse of modified Hessian
    if nargout > 2 % Inverse Hessian vector required
        H_vec = (p_m + p_chi.*k1.*( pi3.^2.*k1 + 2.*eta2.*pi2));
        
        % Improve conditioning of H before inversion
        tau = max(0, sqrt(eps) - (min(H_vec)));
        H_vec = H_vec + tau;  % make Hessian sufficiently positive definite
        Hinv_vec = 1./H_vec;
    end
end
end