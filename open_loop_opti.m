function [j] = open_loop_opti(vehicle, solv)
% Function: optimise the power split between the gas turbine and 
% electric motor of a parallel hybrid propulsion system to meet the drive 
% power demand for a given flight profile. Open loop implementation of the 
% energy management algorithm in the paper "Predictive energy management 
% for hybrid electric aircraft propulsion systems" . 
% 
% Input: 
%        - vehicle: structure of parameters
%        - solv: solver type. 1: CVX, 2: ADMM
% Output: 
%         - j: number of ADMM iterations (return NaN for CVX)
%         - plots of the results for debug purpose

%% Get parameters
m_fuel_0 = vehicle.m_fuel_0;  % initial total fuel mass (kg)
m_0_0 = vehicle.MTOW;         % initial mass of the aircraft (kg)
E_low = vehicle.E_low;        % lower bound on battery energy /shaft (MJ)
E_up = vehicle.E_up;          % uppper bound on battery energy /shaft (MJ)
P_gt_low = vehicle.P_gt_low;  % lower bound on g.t. power /shaft (MW)
P_gt_up = vehicle.P_gt_up;    % upper bound on g.t. power /shaft (MW) 
P_b_low = vehicle.P_b_low;    % lower bound on bat. chem. power /shaft (MW)
P_b_up = vehicle.P_b_up;      % upper bound on bat. chem. power /shaft (MW)
ph_low = vehicle.ph_low;      % lower bound on g.t. fuel rate /shaft (kg/s)
ph_up = vehicle.ph_up;        % upper bound on g.t. fuel rate /shaft (MW)
d = vehicle.d;                % time step (s)
RU = vehicle.RU;              % ratio R_bat/U_bat^2*10^6
n_mot = vehicle.n_mot;        % number of shafts (i.e. propulsion systems)
Phi = vehicle.Phi;            % Lx1 vector of '1' (L is the problem dim.)
t = vehicle.t;                % time vector (s) 
m_drv = vehicle.m_drv;        % drive power coefficients       
n_drv = vehicle.n_drv;
q_drv = vehicle.q_drv;        
alpha = vehicle.alpha;        % e.m. loss map coefficients
beta = vehicle.beta;          % fuel map coefficients
n_iter = vehicle.n_iter;      % max # of ADMM iterations at each time step
sigma1 = vehicle.sigma1;      % ADMM penalty parameter for the 1st constr.      
sigma2 = vehicle.sigma2;      % ADMM penalty parameter for the 2nd constr.
sigma3 = vehicle.sigma3;      % ADMM penalty parameter for the 3rd constr.
sigma4 = vehicle.sigma4;      % ADMM penalty parameter for the 4th constr.
sigma5 = vehicle.sigma5;      % ADMM penalty parameter for the 5th constr.
sigm = vehicle.sigm;          % frequency of sigmas update
tol = vehicle.tol;            % stoping criterion for the residuals 
tol_abs = vehicle.tol_abs;    % e_abs in paper 

%% Define
L = vehicle.L;                % problem dimension (L = length(t))
Psi= d*tril(ones(L));         % lower triangular matrix for dynamic update
E_0_0 = E_up; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CVX Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (solv == 1)
j = NaN;

if (alpha(1, 3) == 0)  % linear em loss map 
 cvx_begin
  cvx_precision best
   variables P_b(L) m(L)
   minimize(sum(n_mot*d*(Psi\(Phi*m_0_0-m)/n_mot-beta(:,1))./beta(:, 2)))
   subject to
    beta(:,2).*((m_drv.*[m_0_0;m(1:end-1)].^2+n_drv.*[m_0_0; m(1:end-1)]...
    + q_drv) - (-RU*P_b.^2 + P_b-alpha(:,1))./alpha(:,2)) + beta(:,1) ...
    - Psi\(Phi*m_0_0-m)/n_mot  <= 0
    P_b_low <= P_b <= P_b_up
    P_gt_low*Phi <= (Psi\(Phi*m_0_0-m)/n_mot...
                                   - beta(:, 1))./beta(:, 2) <= P_gt_up*Phi
    E_low*Phi <= Phi*E_0_0 - cumsum(P_b)*d <= E_up*Phi
 cvx_end
else  % quadratic em loss map
 cvx_begin
  cvx_precision best
   variables P_b(L) m(L)
   minimize(sum(n_mot*d*(Psi\(Phi*m_0_0-m)/n_mot-beta(:,1))./beta(:, 2)))
   subject to
    beta(:,2).*((m_drv.*[m_0_0;m(1:end-1)].^2+n_drv.*[m_0_0;m(1:end-1)]...
    + q_drv) - (-alpha(:, 2) + sqrt(alpha(:, 2).^2 ...
    - 4*alpha(:,3).*(alpha(:,1) - (-RU*P_b.^2 + P_b))))./(2*alpha(:,3)))...
    + beta(:, 1) - Psi\(Phi*m_0_0-m)/n_mot  <= 0
    P_b_low <= P_b <= P_b_up
    P_gt_low*Phi <= (Psi\(Phi*m_0_0-m)/n_mot...
                                   - beta(:, 1))./beta(:, 2) <= P_gt_up*Phi
    E_low*Phi <= Phi*E_0_0 - cumsum(P_b)*d <= E_up*Phi
 cvx_end
end

%% Compute power split
E = Phi*E_0_0 - cumsum(P_b)*d; 
P_eng= (Psi\(Phi*m_0_0-m)/n_mot - beta(:, 1))./beta(:, 2);
P_em = inv_h(inv_f_P_b(P_b, RU), alpha);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADMM Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(solv == 2) 
    
    %% Initialisation (see first set of eq. in p7, second column) 
    P_b = P_b_up;  % warm start (initialise with CDCS solution)
    zeta = P_b; 
    xi = ph_low;  % ph_low = f(P_gt_low*Phi, beta)
    ph = xi; 
    E = trunc(E_0_0*Phi - cumsum(zeta(1:end))*d, E_low, E_up); % SOC update
    m = m_0_0*Phi - cumsum(xi(1:end))*d*n_mot;  % mass update
    pwr = m_drv.*[m_0_0;m(1:end-1)].^2 + n_drv.*[m_0_0;m(1:end-1)] + q_drv;
    chi= trunc(xi - f_phi(pwr, P_b, alpha, beta, RU), 0*Phi, inf*Phi); 
    lambda1=0*Phi;  % multipliers 
    lambda2=lambda1;
    lambda3=lambda1;
    lambda4=lambda1;
    lambda5=lambda1;
    
    %% Precompute commonly used variables
    Psi_t = Psi'; 
    m_0_0_Phi = m_0_0*Phi; 
    E_0_0_Phi = E_0_0*Phi; 
    Inv1 = ((sigma1 + sigma4)*eye(L) + sigma2*Psi_t*Psi*n_mot^2)^(-1); 
    Inv2 = (sigma5*eye(length(Phi)) + sigma3*Psi_t*Psi)^(-1); 
    o = 0*Phi;
    k1_ = 1./alpha(:,2); 
    c=[o ; m_0_0_Phi; o ; E_0_0_Phi; o];
    % Structure of parameters for optimisation
    par.alpha = alpha; par.n_mot = n_mot;  
    par.beta = beta;  par.m_drv = m_drv; par.n_drv = n_drv; 
    par.q_drv = q_drv;  par.d = d; par.RU = RU; par.any = any(alpha(1,3)); 

    %% ADMM iteration update (see first set of eq. in p. 7, first column)
    j = 1;
    res_r = 1e+5; res_s = 1e+5;    % primal and dual residual norm
    eps_r_abs = 0; eps_s_abs = 0;  % tol. on primal and dual res. norm

    while (j<n_iter && (res_r>eps_r_abs+tol || res_s>eps_s_abs+tol))
        
        % Store previous values for dual residual update
        x_ = [chi'; xi'; zeta'; E'; ph'];
        
        % Slack variables update
        pwr = m_drv.*[m_0_0;m(1:end-1)].^2+n_drv.*[m_0_0;m(1:end-1)]+q_drv;
        chi = trunc(xi-f_phi(pwr,P_b,alpha,beta,RU)-lambda1,0*Phi,inf*Phi); 
        
        xi = Inv1*(-Phi*d*n_mot + sigma1*(chi...
        + f_phi(pwr, P_b, alpha, beta, RU) +lambda1) ...
        -sigma2*n_mot*d*cumsum(m(1:end)-m_0_0_Phi(1:end)+lambda2(1:end),...
        'reverse') + sigma4*(ph - lambda4));  
        % Remark 1: d*cumsum(.) replaces multiplication by Psi for faster
        % computation. 
        % Remark 2: 'reverse' keyword is used for Psi^t
        
        zeta = Inv2*(-sigma3*d*cumsum(E(1:end) - E_0_0_Phi(1:end)...
        + lambda3(1:end),'reverse') + sigma5*(P_b - lambda5)); 
        
        % Battery SOC update
        Psi_zeta = [cumsum(zeta(1:end))*d]; 
        E = trunc(E_0_0_Phi -  Psi_zeta -lambda3,E_low, E_up); 

        % Fuel rate update
        ph = trunc(xi+lambda4,ph_low, ph_up);
        
        % Optimise P_b with Newton's method
        par.chi = chi ; par.xi = xi; par.lambda1 = lambda1;
        par.lambda5 = lambda5; par.zeta = zeta; par.sigma1 = sigma1; 
        par.sigma5 = sigma5; par.m_0_0 = m_0_0; par.m = m; 
        par.sum1 = chi - xi + lambda1; par.sum2 = zeta + lambda5; 
        par.P_b_low = P_b_low; par.P_b_up = P_b_up; par.P_drv = pwr; 
        
        P_b = get_optimum_P_b(P_b, par);
        P_b = trunc(P_b, P_b_low, P_b_up);  % P_b_low <= P_b <= P_b_up  
        
        % Optimise m with Newton's method
        par.P_b = P_b; par.sigma2 = sigma2; par.lambda2 = lambda2; 
        Psi_xi = cumsum(xi(1:end))*d*n_mot;  par.g__ = inv_f_P_b(P_b, RU);
        par.sum3 = Psi_xi - m_0_0_Phi + lambda2; 
        par.g_ = inv_h(inv_f_P_b(P_b, RU), alpha);
        
        m = get_optimum_m(m, par);
        
        % Update multipliers
        pwr = m_drv.*[m_0_0;m(1:end-1)].^2+n_drv.*[m_0_0;m(1:end-1)]+q_drv; 
        fphi = f_phi(pwr, P_b, alpha, beta, RU);
        lambda1 = lambda1 + chi -xi + fphi;
        lambda2 = lambda2 + m - m_0_0_Phi + Psi_xi ;
        lambda3 = lambda3 + E - E_0_0_Phi + Psi_zeta;
        lambda4 = lambda4 + xi - ph ;
        lambda5 = lambda5 + zeta - P_b;
        
        %% Residuals
        
        % Primal: b(u) + Bx = c
        b_u = [fphi; m; -P_b; o; o];
        B_x = [chi - xi; n_mot*d*cumsum(xi(1:end));... 
               zeta; d*cumsum(zeta(1:end)) + E ; xi - ph]; 
        res_scale = max(norm(b_u), max(norm(B_x), norm(c))); 
        res_rel = (b_u + B_x -c)/res_scale; 
        eps_r_abs = tol_abs*sqrt(5*L)/res_scale;
        res_r = norm(res_rel);
        
        % Dual : db'*Rs*B*(x_ - x)
        x = [chi'; xi'; zeta'; E'; ph']; 
        dx = x_ - x;
        gradf_m = beta(:, 2).*(2*m_drv.*[0; m(1:end-1)]...
                  + n_drv.*[0; ones(size(m(1:end-1)))]); 
        
        if (par.any)  % quadratic em loss map
            gradf_P_b = beta(:, 2).*k1_.*(2*par.RU*P_b - 1);
        else          % linear loss map
            gradf_P_b = beta(:,2).*(2*par.RU*P_b-1)./sqrt(alpha(:,2).^2 ...
                - 4*alpha(:,3).*(alpha(:,1) - (-par.RU*P_b.^2 + P_b))); 
        end 

        s_scale = norm([gradf_m.*lambda1 + lambda2;...
                        gradf_P_b.*lambda1 - lambda3]); 
        s_rel = [sigma1*gradf_m.*dx(1,:)' - sigma1*gradf_m.*dx(2,:)'...
        +sigma2*n_mot*d*cumsum(dx(2,1:end)'); sigma1*gradf_P_b.*dx(1,:)'...
        -sigma1*gradf_P_b.*dx(2,:)' - sigma3*dx(3,:)']/s_scale;   
        res_s = norm(s_rel);
        eps_s_abs = tol_abs*sqrt(2*L)/s_scale;
    
        %% Recompute sigmas 
        res_up = max(norm(res_rel*res_scale), norm(s_rel*s_scale)); 
        if (mod(j,sigm) == 0 && 10 < res_up)
            
            % Recompute sigmas
            [sigma1, sigma2, sigma3, sigma4, sigma5, mult] = sig(sigma1,...
                           sigma2, sigma3, sigma4, sigma5, res_rel, s_rel);
            
            % Rescale multipliers
            lambda1 = lambda1*mult(1);    
            lambda2 = lambda2*mult(2);
            lambda3 = lambda3*mult(3);
            lambda4 = lambda4*mult(4);
            lambda5 = lambda5*mult(5);
    
            Inv1 = ((sigma1+sigma4)*eye(L)+sigma2*Psi_t*Psi*n_mot^2)^(-1); 
            Inv2 = (sigma5*eye(L) + sigma3*Psi_t*Psi)^(-1);
        end 
        j = j+1; 
    end 
    
    %% Compute power split 
    P_eng= (ph - beta(:, 1))./beta(:, 2);     % gt power /shaft (MW)
    P_em = inv_h(inv_f_P_b(P_b, RU), alpha);  % battery power /shaft (MW)

end 
 
%% Hybridisation ratio
P_drv = m_drv.*[m_0_0; m(1:end-1)].^2 + n_drv.*[m_0_0; m(1:end-1)] + q_drv; 
H_ratio = P_em./P_drv;

%% Plots 
figure()
plot(t, H_ratio)
ylabel('Hybridisation ratio P_{em}/P_{drv} (-)')
xlabel('Time (s)')

figure()
plot(t, P_em + P_eng)
hold on
plot(t, P_drv)
legend('P_{em} + P_{gt}', 'P_{drv}')
ylabel('Power Check (MW)')
xlabel('Time (s)')

figure()
plot(t, [m_0_0; m(1:end-1)] - m_0_0 + m_fuel_0)
ylabel('Fuel tank level (kg)')
xlabel('Time (s)')
hold on

figure()
plot(t, [E_0_0; E(1:end-1)])
hold on
plot(t, ones(size(t))*E_low)
ylabel('SOC (MJ)')
xlabel('Time (s)')

figure()
plot(t, P_em)
hold on
plot(t, P_eng)
plot(t, P_drv)
legend('P_{em}', 'P_{gt}', 'P_{drv}')
ylabel('Power profiles (MW)')
xlabel('Time (s)')
end 