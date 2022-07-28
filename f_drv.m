function [eta2, eta1, eta0] = f_drv(v, v_dot, gamma, gamma_dot,...
                                                  Cd, Cl, g, S, rho, nprop)
% Function: Compute drive power coefficients in the expression relating
% drive power P_drv to the aircraft mass m and flight parameters:
%                   P_drv = eta2.*m.^2 + eta1.*m + eta0
% (see equation 26 in the paper "Predictive energy management for hybrid
% electric aircraft propulsion systems"). 
% Note that P_drv in this form is the drive power per shaft in MW. 
% 
% Input:
%        - v: predicted velocity profile 
%        - v_dot: predicted acceleration profile
%        - gamma: predicted flight path angle profile
%        - gamma_dot: predicted rate of flight path angle profile
%        - Cd: vector of drag coefficients
%        - Cl: vector of lift coefficients
%        - g: gravity acceleration
%        - S: characteristic surface
%        - rho: air density
%        - nprop: number of shafts (i.e. number of gt + em systems)
%
% Output: eta2, eta1, eta0: power coefficients appearing in 
% P_drv = eta2.*m.^2 + eta1.*m + eta0

% Get drag and lift coefficients
b0 = Cl(1); b1 = Cl(2);
a0 = Cd(1); a1 = Cd(2); a2=Cd(3); 

% Eta coefficients as in eq. 15 of paper
eta2 = (2*a2*(v.*gamma_dot + g*cos(gamma)).^2)./(b1^2*rho*S.*v);
eta1 = v_dot.*v + g*sin(gamma).*v - 2*a2*b0*v.*(v.*gamma_dot ...
    + g*cos(gamma))/b1^2 + a1/b1*(v.*gamma_dot + g*cos(gamma)).*v; 
eta0 = 0.5*rho*S.*v.^3*(a0 + a2*b0^2/b1^2 - a1*b0/b1);

% Conversion from total W -> MW per shaft
eta2 = eta2*1e-6/nprop;
eta1 = eta1*1e-6/nprop;
eta0 = eta0*1e-6/nprop; 

end