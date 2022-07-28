function beta = interp_gt(v, h, n)
% Function: Interpolate the gas turbine fuel map coefficients from gt data. 
% The fuel map is given by:
%
%               m_fuel = beta1*P_em + beta0
%
% The gas turbine data is organised in 36 blocks of 11 data points for 
% which the height (altitudem) and mach number (MachNo) are onstant. For 
% each block, fuel map coefficients are fitted to the fuel consumption 
% (Wffkgs) and shaft power (ShaftPowerW) entries (11 data points  /block).
% This allows to generate a grid of fuel coefficients associated to each 
% pair of height and mach no. 
%
% The coeff. are then interpolated for the input estimated 
% velocity and height profile from the data grid.
%
% Input: 
%        - v: estimated velocity profile
%        - h: estimated height profile
%
% Ouput: 
%        - beta: vector of fuel map coeff. beta0 (kg/s), beta1 (kg/s/MW)


%% Load data
load('allison.mat')
altitude = Allison250PerformanceDeckFullRange.altitudem;  % altitude (m)
mach = Allison250PerformanceDeckFullRange.MachNo;  % Mach number (-)
power = Allison250PerformanceDeckFullRange.ShaftPowerW;  % power output (W)
fuel = Allison250PerformanceDeckFullRange.Wffkgs;  % fuel consumpt° (kg/s)
 
%% Fit loss map coefficients for a given alt and Mach #
N = length(v); 

offset = 11;  % number of data points per block
N_blck = length(altitude)/offset;  % numnber of blocks
Beta = zeros(4, 9, n+1);
Altitude = zeros(4, 9);
Mach = zeros(4, 9);
k = 0; % initialise block number

% Explore the data block by block
for i=1:4
    for j=1:9
        k = k + 1;

        % Assemble meshgrid
        Altitude(i, j) = altitude(1+offset*(k-1)); 
        Mach(i, j) = mach(1+offset*(k-1));
        
        % Collect data for the current block
        P_gt = power(1+offset*(k-1):offset*k)/10^6;  % convert in MW
        m_fuel = fuel(1+offset*(k-1):offset*k);
        
        % Least square fit
        x= [ones(offset, 1) P_gt]\m_fuel;  % first order polynomial fit
        Beta(i, j, :) = x';
    end 
end 

%% Mach No for estimated profile
load('Sound_Data.mat')
c = interp1(sound_data.alt, sound_data.c, h); 
M = v./c; 

%% Interpolate coefficients for estimated profile 
beta = zeros(N, n+1);

for i=1:n+1
     % Interp.: 'spline' method allows extrapolation beyond data range
    beta(:, i) = interp2(Altitude, Mach, Beta(:, :, i), h, M, 'makima');  
end 
end