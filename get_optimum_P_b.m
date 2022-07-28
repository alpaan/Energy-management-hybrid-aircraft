function P_b = get_optimum_P_b(P_b, par)
% Function:  Solve the optimisation problem in the P_b update step of the 
% ADMM algorithm using Newton's method with a modified Hessian and
% backtracking line search (Wolfe-Powell step-size rule). 
% Input: 
%        - P_b: current estimate of battery chemical power per shaft (MW)
%        - par: structure containing parameters and variables
%
% Output: P_b: updated battery chemical power per shaft (MW)

% Handle to cost function, gradient and Hessian
cost_grad_hess_P_b = @(x) f_eval_P_b(x, par.sum1, par.P_drv, par.sum2, par);

% Projected Line Search P_b
tol = 1e-7;
for j = 1:1
    [P_b, flag] = newton_line_search(P_b, cost_grad_hess_P_b, tol);
     P_b = trunc(P_b, par.P_b_low, par.P_b_up);  % P_b_low <= P_b <= P_b_up
    if flag == -1
        P_b = par.P_b_low;
    elseif flag == 1
        break;
    end
end
end