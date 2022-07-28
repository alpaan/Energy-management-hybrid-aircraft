function [m_dot] = f(P_gt, beta)
% Function: Compute the rate of fuel mass (kg/s) from the engine power (MW)
% Input:
%          - P_gt: gas turbine power (MW)
%          - beta: vector of fuel map coefficients b0, b1 (-, MW^-1)
% Output:
%          - m_dot: fuel mass (kg/s)        

m_dot = beta(:, 2).*P_gt + beta(:, 1);  
end

