function m = get_optimum_m(m, par)
% Function:  Solve the optimisation problem in the m update step of the 
% ADMM algorithm using Newton's method with a modified Hessian and
% backtracking line search (Wolfe-Powell step-size rule). 
% Input: 
%        - m: current estimate of total aircraft mass (kg)
%        - par: structure containing parameters and variables
%
% Output: m: updated total aircraft mass (kg)


cost_grad_hess_m = @(x) f_eval_m(x, par.sum1, -par.g_, par.sum3, par);

% Line Search over m
tol = 1e-6;
for j = 1:1
    [m, flag] = newton_line_search(m, cost_grad_hess_m, tol);
    if flag == 1
        break;
    end 
end

end