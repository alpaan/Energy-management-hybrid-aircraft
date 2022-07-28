%% Energy management for hybrid electric aircraft
% Optimise the power split between the gas turbine and 
% electric motor of a parallel hybrid propulsion system to meet the drive 
% power demand for a given flight profile. Open loop implementation of the 
% energy management algorithm in the paper "Predictive energy management 
% for hybrid electric aircraft propulsion systems".
%
% Creation: Martin Doff-Sotta (October 2019) 
% Last update: Martin Doff-Sotta (July 2022)
%
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solver parameters
solve_type = 2;             % solver (1: CVX, 2: ADMM, 3: CDCS)
d = 60;                     % sample time (s)
T = 3600;                   % duration of flight (s) 
n_iter = 20000;             % max # of ADMM iterations at each time step  
sigma1 = 5e+1;              % ADMM penalty parameter for the 1st constr.             
sigma2 = 3.6905e-1;         % ADMM penalty parameter for the 2nd constr.
sigma3 = 6.9621e-5;         % ADMM penalty parameter for the 3rd constr.
sigma4 = 2.02857e-1;        % ADMM penalty parameter for the 4th constr.
sigma5 = 8.29e-1;           % ADMM penalty parameter for the 5th constr.
sigm = 50;                  % frequency of sigmas update 
                            % (put sigm > n_iter to keep sigmas fixed) 
tol = 5e-8;                 % stoping criterion for the residuals 
tol_abs = 0;                % e_abs in paper 

% Physical parameters
n_mot = 4;                  % number of propulsion systems 
m_bat = 8000;               % battery mass for all systems (kg)
dens = 0.875;               % battery density (MJ/kg)
SOC_max = dens*m_bat/n_mot; % max battery stored energy per system (MJ)
R_bat = 3.5e-2;             % battery eqvl circuit resistance (ohm)
U_bat = 1500;               % battery eqvl circuit open circuit voltage (V)
m_fuel_0 = 4000;            % initial fuel mass (kg)
E_low = 0.2*SOC_max;        % lower bound on battery energy per system (MJ)
E_up = 0.85*SOC_max;        % upper bound on battery energy per system (MJ)
P_em_low = 0;               % lower bound on e.m. power per shaft (MW)
P_em_up = 5;                % upper bound on e.m. power per shaft (MW) 
P_gt_low = 0;               % lower bound on gas turbine power (MW)
P_gt_up = 5;                % upper bound on gas turbine power (MW)
compCl = [0.11 0.43];       % lift coefficients: [b0(-); b1(deg^-1)]
compCd = [0.0282 0.0042...
              5.3613e-04];  % drag coeff.: [a0(-); a1(deg^-1); a2(deg^-2)]
initTAS = 131.9263;         % initial true airspeed (m/s)
g = 9.81;                   % acceleration due to gravity (m/s^2)
S = 77.3;                   % characteristic surface area (m^2)
M = 42000;                  % aircraft initial mass (kg)
eta_w = 0.5;                % battery charge up efficiency (-) 
w_mil = 0;                  % 0: normal mode, 1: windmilling (ADMM only)
if (w_mil == 1 && (solve_type == 1 || solve_type == 3))
    disp('Solver error: windmilling works only with ADMM.')
    disp(' Try setting solve_type = 2 or w_mil = 0.')
    return;
elseif (w_mil == 1 && solve_type == 2)
    P_em_low = - P_em_up; 
end 

% Load gas turbine data
load("allison_speed.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Trajectory planning %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flight parameters 
vtas_up = 190;                           % cruise TAS (m/s)
h_up = 12750;                            % cruise altitude (m)
rho_air = get_rho(h_up);                 % air density at ref alt. (kg/m^3) 

% Flight profile generation
vp=[initTAS  vtas_up...
    vtas_up initTAS]';                   % TAS data points (m/s)
tvp = [0 500 2000 T]';                   % TAS time data points (s)
hp = [0 h_up h_up 0];                    % height data points (m)
thp = [0 1000 2000 T]';                  % height time data points (s)
t=(0:d:T)';                              % time at EM sampling time (s)
v = interp1(tvp, vp, t, 'pchip');        % interp. of TAS profile (m/s)
height = interp1(thp, hp, t, 'pchip');   % interp. of height profile (m)
h_dot=diff(height)./diff(t);             % vertical speed (m/s) 
gamma=asin(h_dot./v(1:end-1));
gamma=[gamma; gamma(end)];               % flight path angle (rad) (level)
v_dot = diff(v)./diff(t); 
v_dot = [v_dot; v_dot(end)];             % acceleration profile (m/s2)
gamma_dot = diff(gamma)./diff(t); 
gamma_dot = [gamma_dot; gamma_dot(end)]; % rate of flight path angle (m/s2)
height = h_up*ones(size(v));             % height profile (m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Power profile & losses %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(t); 
Phi = ones(length(t), 1);

% Get drive power coeff. in P_drv = m_drv.*m.^2 + n_drv.*m + q_drv
[m_drv, n_drv, q_drv] = f_drv(v, v_dot, gamma, gamma_dot, compCd, compCl,...
                                    g, S, rho_air, n_mot); % (MW per shaft)

% Interpolate gt fuel coeff. in ph = beta1*Peng + beta0 (1st order fit)
%beta = interp_gt(v, height, 1);
beta = [0.0327 0.0821].*ones(N, 2);

% Estimated drive power when gas turbine is iddle (MW per shaft)
m_est = M - [0;d*cumsum(n_mot*f(P_gt_low*Phi(1:end-1),beta(1:end-1)))];
P_drv_est = m_drv.*m_est.^2 + n_drv.*m_est + q_drv;

% Interp. e.m. loss map coeff. in Pc = alpha2*P_gt^2 + alpha1*P_gt + alpha0 
w_gt_est = interp3(gas_turbine.h, gas_turbine.v, gas_turbine.P_gt, ...
gas_turbine.w, height, v, P_drv_est, 'spline'); % Shaft rotation speed 
%alpha = interp_em(w_gt_est, 1);  % order of fit: 1 or 2. 
alpha = [0 1 0].*ones(N, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Power constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redefine electric motor output power (mechanical) to ensure that battery 
% loss map f_P_b(.) is real-valued 
P_b_up = zeros(N, 1);
P_b_low = zeros(N, 1);
P_em_up = P_em_up*Phi;
P_em_low = P_em_low*Phi;
RU = R_bat/U_bat^2*10^6; 
for i=1:length(t)
    % Upper bound on battery chemical power P_b
    x_P_em=roots([-4*RU*alpha(i, 2)  1-4*RU*alpha(i, 1)]); 
    P_em_up(i,1) = min(P_em_up(i,1), max(x_P_em));
    P_b_up(i,1)= f_P_b(h(P_em_up(i,1), alpha(i,:)), RU); 
    
    % Lower bound on battery chemical power P_b
    if (any(alpha(1,3)))
        P_em_low(i, 1) = max(P_em_low(i,1), -alpha(i, 2)./(2*alpha(i, 3)));
        P_b_low(i,1)= f_P_b(h(P_em_low(i,1), alpha(i,:)), RU);
    else
        P_b_low(i,1) = f_P_b(h(P_em_low(i,1), alpha(i,:)), RU);
    end 
end

% Bounds on electric motor input (electrical) power P_c
P_c_low = h(P_em_low, alpha);
P_c_up = h(P_em_up, alpha);

% Bounds on fuel consumption
ph_low = f(P_gt_low*Phi, beta);
ph_up = f(P_gt_up*Phi, beta);

% Bounds on commanded drive power
P_drv_low = max(P_em_low + P_gt_low);
P_drv_up = min(P_em_up + P_gt_up);

% Windmilling 
if (w_mil == 1)
    ind = find(P_drv_est < 0); 
    P_b_low(ind,1)= f_P_b(h(P_drv_est(ind,1)*eta_w, alpha(ind,:)), RU); 
    P_b_up(ind,1)= f_P_b(h(P_drv_est(ind,1)*eta_w, alpha(ind,:)), RU);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Open loop optimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi= d*tril(ones(N));         % lower triangular matrix for dynamic update
E_0_0 = E_up; 
m_0_0 = M; 

% CVX solver 
if (solve_type == 1)

if (alpha(1, 3) == 0)  % linear em loss map 
 cvx_begin
  cvx_precision best
   variables P_b(N) m(N)
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
   variables P_b(N) m(N)
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

% Compute power split
E = Phi*E_0_0 - cumsum(P_b)*d; 
P_eng= (Psi\(Phi*m_0_0-m)/n_mot - beta(:, 1))./beta(:, 2);
P_em = inv_h(inv_f_P_b(P_b, RU), alpha);


% ADMM solver
elseif(solve_type == 2) 
    
    % Initialisation (see first set of eq. in p7, second column) 
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
    
    % Precompute commonly used variables
    Psi_t = Psi'; 
    m_0_0_Phi = m_0_0*Phi; 
    E_0_0_Phi = E_0_0*Phi; 
    Inv1 = ((sigma1 + sigma4)*eye(N) + sigma2*Psi_t*Psi*n_mot^2)^(-1); 
    Inv2 = (sigma5*eye(length(Phi)) + sigma3*Psi_t*Psi)^(-1); 
    o = 0*Phi;
    k1_ = 1./alpha(:,2); 
    c=[o ; m_0_0_Phi; o ; E_0_0_Phi; o];

    % Structure of parameters for optimisation
    par.alpha = alpha; par.n_mot = n_mot;  
    par.beta = beta;  par.m_drv = m_drv; par.n_drv = n_drv; 
    par.q_drv = q_drv;  par.d = d; par.RU = RU; par.any = any(alpha(1,3)); 

    % ADMM iteration update (see first set of eq. in p. 7, first column)
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
        
        % Primal residual: b(u) + Bx = c
        b_u = [fphi; m; -P_b; o; o];
        B_x = [chi - xi; n_mot*d*cumsum(xi(1:end));... 
               zeta; d*cumsum(zeta(1:end)) + E ; xi - ph]; 
        res_scale = max(norm(b_u), max(norm(B_x), norm(c))); 
        res_rel = (b_u + B_x -c)/res_scale; 
        eps_r_abs = tol_abs*sqrt(5*N)/res_scale;
        res_r = norm(res_rel);
        
        % Dual residual: db'*Rs*B*(x_ - x)
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
        eps_s_abs = tol_abs*sqrt(2*N)/s_scale;
    
        % Recompute penalty parameters
        res_up = max(norm(res_rel*res_scale), norm(s_rel*s_scale)); 
        if (mod(j,sigm) == 0 && 10 < res_up)
            
            % Update sigmas
            [sigma1, sigma2, sigma3, sigma4, sigma5, mult] = sig(sigma1,...
                           sigma2, sigma3, sigma4, sigma5, res_rel, s_rel);
            
            % Rescale multipliers
            lambda1 = lambda1*mult(1);    
            lambda2 = lambda2*mult(2);
            lambda3 = lambda3*mult(3);
            lambda4 = lambda4*mult(4);
            lambda5 = lambda5*mult(5);
    
            Inv1 = ((sigma1+sigma4)*eye(N)+sigma2*Psi_t*Psi*n_mot^2)^(-1); 
            Inv2 = (sigma5*eye(N) + sigma3*Psi_t*Psi)^(-1);
        end 
        j = j+1; 
    end 
    ADMM_iterations = j
    
    % Compute power split 
    P_eng= (ph - beta(:, 1))./beta(:, 2);     % gt power /shaft (MW)
    P_em = inv_h(inv_f_P_b(P_b, RU), alpha);  % battery power /shaft (MW)

% CDCS solver 
elseif(solve_type == 3)
    
    % Initialisation 
    P_drv = zeros(size(t));
    P_eng = zeros(size(t));
    P_em = zeros(size(t));
    m = zeros(size(t));
    E = zeros(size(t));
    m(1) = M;
    E(1) = E_up;
    
    % At each time step, battery is used at its max capacity
    for k=1:length(t)
        P_drv(k) = m_drv(k).*m(k).^2+n_drv(k).*m(k)+q_drv(k); % drive power 
        if (E(k) > E_low)  % if battery is not depleted 
            if (P_drv(k) >= P_em_up)  % if power demand saturates em power
                P_em(k) = P_em_up;  % em operates at max capacity 
            else
                P_em(k) = P_drv(k);
            end
            
            % Correct excess battery depeletion due to finite time step
            if(E(k) - d*f_P_b(h(P_em(k), alpha(k,:)), RU) < E_low)
                P_em(k) = inv_h(inv_f_P_b((E(k)-E_low)/d, RU), alpha(k,:));
            end 
            P_eng(k) = P_drv(k) - P_em(k);

        else 
            if (P_drv(k) > 0)
                P_eng(k) = P_drv(k);
                P_em(k) = 0;
            else
                P_eng(k) = 0;
                P_em(k) = 0;
            end
        end 
        
        % Mass update
        m_dot = f(P_eng(k), beta(k, :))*n_mot;  % fuel rate 
        m(k+1) = m(k) - d*m_dot;  % whole aircraft mass update
        
        % Battery SOC update (per battery)
        E(k+1) = E(k) -d*f_P_b(h(P_em(k), alpha(k,:)), RU);
    end
    
    % Shift results to match ADMM and CVX outputs 
    E = E(2:end); 
    m = m(2:end);
    
end 
 
% Hybridisation ratio
P_drv = m_drv.*[m_0_0; m(1:end-1)].^2 + n_drv.*[m_0_0; m(1:end-1)] + q_drv; 
H_ratio = P_em./P_drv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burnt_fuel = M - m(end)

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