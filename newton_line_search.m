function [xi, flag] = newton_line_search(x0, fun_grad_hess, tol)
% Function: Newton backtracking line-search method (Wolfe-Powell rule)
% Input:  
%         - x0: initial guess
%         - fun_grad_hess: function computing cost, gradient and Hessian for
%         the problem at hand
%         - tol: algorithm tolerance
% Output:
%         - xi: optimal value
%         - flag: indicate if algorithm has converged (1) or not (-1)

% Wolfe-Powel update rule parameters (see Nocedal for recommended value)
c1 = 1e-4;
c2 = 0.9;

% Get cost, gradient and (inverse modified) Hessian
[f0, g0, Hinv0] = fun_grad_hess(x0);
if ~isreal(f0)
    flag = -1; 
    xi = NaN;
    return
end

% Search Direction (could include regularisation in the hessian)
p = -Hinv0.*g0; steep_desc = false;
slope0 = g0'*p;
if slope0 >= 0	% Slope is uphill
    p = -g0;	% Use steepest descent instead
    slope0 = g0'*p;
    steep_desc = true;
end

% Performs search for point satisfying the two strong Wolfe conditions:
%	f(x + a * p) <= f(x) + c1 a p'.g(0)		    (WC1)
%	|p' . g(x + a * p)| >= c2|p' . g'(x)|		(WC2b)
i = 1;
a_max = 5;
ai = 0; ai_1 = 0; fi_1 = 0; slopei_1 = 0; 
while 1 
    if i == 1
        if ~steep_desc
            % Newton method, use ai = 1 as first trial
            ai = 1;
        else
            ai = 0;
        end
        ai_1 = 0;
        fi_1 = f0;
        slopei_1 = slope0;
    end
    % Evaluate phi(a_i)
    xi = x0 + ai.*p;
    [fi, gi] = fun_grad_hess(xi);
    % Test Wolfe Condition 1
    WC1 = fi <= f0 + c1.*ai.*slope0;
    fun_incr = fi >= fi_1 && i ~= 1;
    if (~WC1) || fun_incr
        % WC1 doesn't hold or function increase, enter zoom
        [xi, gi] = zoom(ai_1, ai, fi_1, fi, slopei_1);
        break;
    end
    % WC1 holds, test WC2
    slopei = gi'*p;
    %     WC2a = slopei >= c2.*slope0;
    WC2b = abs(slopei) <= c2.*abs(slope0);
    if WC2b
        break;
    end
    if slopei >= 0
        [xi, gi] = zoom(ai, ai_1, fi, fi_1, slopei);
        break;
    end
    
    % Update ai and fi
    ai_1 = ai;
    ai = new_a_mono_inc(ai, a_max);
    fi_1 = fi;
    % Update iteration counter
    i = i + 1;
end

% Check for convergence
flag = double(all(abs(gi) < tol .* ones(size(x0))));


% Functions
function [xj, gj] = zoom(a_lo, a_hi, f_lo, f_hi, slope_lo)
    max_zoom_iterations = 2;
    j = 0;

    while 1
        % quadratic interpolant of f values on a_lo, a_hi
        % and find min of interpolant within interval [a_lo, a_hi]
        da = a_hi - a_lo;
        q1 = slope_lo;
        q2 = (f_hi - f_lo - da*slope_lo)/da^2;
        if (q2 <= 0 || ~isreal(q2))
            % Use bisection
            aj = a_lo + 0.5*da;
        else
            % Use min of quadratic
            aj = a_lo - 0.5*q1/q2;
        end
            
        % Evaluate at new a
        xj = x0 + aj*p;
        [fj, gj] = fun_grad_hess(xj);
            
        % Check decrease condition
        WC1 = fj <= f0 + c1.*aj.*slope0;
        fun_incr = fj >= f_lo;
        slopej = gj'*p;
        if (~WC1) || fun_incr
            a_hi = aj;
            f_hi = fj;
        else
            WC2b = abs(slopej) <= c2.*abs(slope0);
            if WC2b
                return;
            end
            if ( slopej*da >= 0 )
                a_hi = a_lo;
                f_hi = f_lo;
            end
            a_lo = aj;
            f_lo = fj;
            slope_lo = slopej;
        end

        % Update iteration count
        j = j + 1;
        if j > max_zoom_iterations
            break;
        end
    end
end
end

function a = new_a_mono_inc(ai, a_max)
    a = min(a_max, 2*ai);
end