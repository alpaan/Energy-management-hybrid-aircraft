function [x_c] = trunc(x, x_low, x_up)
% Function: Projector function: truncates the input x on the interval
% [x_low ; x_up]
% Input: - x: variable to project
%        - x_low, x_up: lower and upper bounds of the interval
% Output: x_c: truncated variable

x_c = min(x_up, max(x, x_low)); 
end

