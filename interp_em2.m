function alpha = interp_em(w_est, n)
% Function: Interpolate the electric motor loss coefficients from data. For
% a second order fit, the loss map is given by:
%
%               P_c = alpha2*P_em^2 + alpha1*P_em + alpha0
%
% For a first order fit, the loss map is: 
%
%               P_c = alpha1*P_em + alpha0, alpha0 >= 0 
%
% Loss map coefficients are fitted to the data at each velocity data point
% solving a least square problem for the 
% second order fit (resp. constrained least square for the first order fit)
% The alpha coefficients are then interpolated for the input estimated 
% velocity vector w_est. with a spline interpolation method.
%
%
% Input: 
%        - w_est: predicted shaft rotation speed (rad/s)
%        - n: order of fitting polynomial (only 1 or 2) 
%
% Ouput: alpha: vector of loss map coefficients alpha0 (MW), alpha1 (-),
% alpha2 (1/MW)


%% Load data
load('Electric_Motor.mat')
mech_power = electric_motor.P_em;  % output (mechanical) power per shaft (MW)
elec_power = electric_motor.P_c;  % input (electrical) power per shaft (MW)
w = electric_motor.w;  % shaft rotation speed (rad/s)
 
%% Fit loss map coefficients to P_em and P_c for a given speed w
N = length(w); 
A = zeros(N, n+1);

for i=1:N  % loop over speed data

    % Prepare data
    P_em = mech_power(:, i);  % P_em vector for i_th speed datapoint
    P_c = elec_power(:, i);  % P_c vector for i_th speed datapoint
    %index = find(P_em < P_em_up);  % keep only data below P_em_up
    %P_em = P_em(index); 
    %P_c = P_c(index);
    
    % Least square
    Y = P_c;  % target vector
    if (n == 1)
        X = [ones(size(P_em)) P_em ];  % fit 1st order poly
        beta = X\Y; 
    elseif (n == 2)
        X = [ones(size(P_em)) P_em P_em.^2];  % fit 2nd order poly
        beta = X\Y;  % least square fit (solves normal equation)
    end 
    
    A(i, :) = beta';  % store loss map coefficients at i_th speed
end 

%% Interpolate coefficients for estimated shaft rotation speed profile 
alpha = zeros(length(w_est), 3);

for i=1:n+1
    alpha(:, i) = interp1(w, A(:,i), w_est, 'spline');  % interpolation: 'spline' method allows extrapolation beyond data range
end 

end