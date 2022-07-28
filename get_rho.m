function rho = get_rho(h)
% Function: compute air density as a function of altitude h
% Input: h: altitude (m)
% Output: rho: air density (kg/m^3)

% International Standard Atmosphere,
p0 = 101325;    % sea level standard atmospheric pressure (Pa)
T0 = 288.15;    % sea level standard temperature (K)
g = 9.80665;    % earth-surface gravitational acceleration (m/s2)
L = 0.0065;     % temperature lapse rate (K/m)
R = 8.31446;    % ideal (universal) gas constant (J/(mol·K))
M = 0.0289652;  % molar mass of dry air (kg/mol)

% Tropospheric model (h < 10km)
rho = p0*M/(R*T0)*(1 - L*h/T0).^(g*M/(R*L) - 1); 
end 