%% Clear workspace
clear
clearvars -except files nn
close all
clc

% Select execution time profiler
profile off
warning('OFF')

%--------------------------------------------------------------------------
% Plot settings
%--------------------------------------------------------------------------

ThreeD      = false;            % Plot 3D results: true / false
Times_3D    = [5, 15, 20, 30];  % Times for which to generate 3D plots [s]
Arr_sum    = false;            % Calculate arriving particles:  true / false
Plot_1D     = false;            % Plot 1D results of middle of the plume: true / false
Plot_BG     = false;            % Show bg in plot: true / false

%--------------------------------------------------------------------------
% Simulation setting
%--------------------------------------------------------------------------

% Restrictions (Check values)
enableOxidation = false;    %?% Enable particle oxidation (true = no oxidization)
enableDioxides  = false;    % Let particles form dioxides: true / false
nMinPerBin = 1;             % Minimal number of particles per bin
nMinPerAngle = 1E7;         % Minimal N per angle (angles with lower N are not computed)
nCollisionsMax = 7;         % Maximum number of collisions per particle

% Smoothing parameters
Smooth_Fraction = 1E4;  % Fraction to smooth % ???
Smooth_Value    = 15;   % Moving average smooth value (default: 15)

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

%?%n_A     = find(M_Array == M_A) - 1;     % place of A atom in the array
%?%n_B     = find(M_Array == M_B) - 1;     % place of B atom in the array
%?%n_O     = find(M_Array == M_O) - 1;     % place of oxygen atom in the array
%?%n_AO    = find(M_Array == M_A+M_O) - 1; % place of AO molecule in the array
%?%n_BO    = find(M_Array == M_B+M_O) - 1; % place of BO molecule in the array
%?%n_BO2   = find(M_Array == M_B+M_O) - 1; % place of BO2 molecule in the array

%?%n_position = [n_A n_B n_O n_AO n_BO n_BO2]; %array of the position in m_p array

% Target
%?%Density_Target = 1;  % Target density

% Energies
E_TiO = 0.117 * CONSTANT.EV;                 % Activation energy for forming TiO [J]
E_O_low = 0.5 * CONSTANT.EV;                 % Lower bound oxidation energy [J]
E_O_high_Sr = 4.418 * CONSTANT.EV; %* Oxidize; % Upper bound oxidation energy Sr [J]
E_O_high_Ti = 2.276 * CONSTANT.EV; %* Oxidize; % Upper bound oxidation energy Ti [J]

%?%t_p = 4;     %number of start particles +1
%?%% t_o = 4;   %number of oxidation states
%?%% t_oo = 2;  %placement oxide particle.
%?%% t_pp = 8;  %total number of particles

%--------------------------------------------------------------------------
% Deposition settings
%--------------------------------------------------------------------------

% Background gas
P_BG_PASCAL = 100 * bgPressure; % Convert pressure in Pa (1 mbar = 10^2 Pa)

% System dimensions
X_Ts = 0.05;    % Target-substrate distance [m]

%--------------------------------------------------------------------------
% Initial plume settings
%--------------------------------------------------------------------------

dis_NPoints =  100; % Number of datapoints in the initial velocity distribution

%--------------------------------------------------------------------------
% Initial distribution
%--------------------------------------------------------------------------

% N_UC = (SPOT_VOLUME * UC_N_TOTAL) / UC_Vol; % * Density_Target)); % Total number of ablated atoms per UC
%?%N_TOTAL = (UC_N_Atoms * N_UC) / UC_N_TOTAL;  % Total number of ablated atoms
