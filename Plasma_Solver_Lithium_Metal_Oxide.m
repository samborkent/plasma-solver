% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plasma Solver v5.0.0 16-02-21
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Author:     Tom Wijnands
%     Modified for two metal oxides by:     Bart Boonstra
% Completely reworked and modified for
%              lithium metal oxides by:     Sam Borkent
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Description:
% 
% Based on 'Numerical modelling of the plasma plume propagation and
% oxidation during pulsed laser deposition of complex oxide thin films', by
% 2020, by T. Wijnands, E.P. Houwman, G. Koster, G. Rijnders, M. Huijben.
% 
% The original script by Tom Wijnads modeled the propagation of an PLD
% plasma plume in 3D as result of ablation of a single crystal TiO2 target,
% initial ablation is not included. Modified by Bart Boonstra to improve
% performance and include two metals (SrTiO3). Modified and refactored by
% Sam Borkent to include two metals with highly different masses (Li2TiO3).
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------------
% Change log:
% --------------------------------------------------------------------------
% 
% 17-02-21:
%   * Edited to include Git support
%   * Changed naming conventions to more closely resemble the names used
%       in the paper and to follow MATLAB style guidelines
%   * Changed maximum velocity to include limits of Li
%   * Wrote the initial velocity distibution as a function
% 
% 16-02-21:
%   * Work on initial velocity distribution
% 
% 12-02-21:
%   * Started refactoring main loop
%   * Replaced looping through angles to index looping
%   * Refactored collision cross section
% 
% 08-02-21:
%   * Starting top-to-bottom to refactor the code to increase readability
%   * Adjust naming conventions and replace Dutch variable names
%   * Change 1/0 bools into true/false bools
%   * Denoting lines with unclear purpose with %?%
%   * Reorder variables and place them in more clear groups
%   * Added distiction between lower/Camel/UPPER case for
%       variables/Parameter/CONSTANTS
%   * Made make dir conditional so it is skipped if dir already exists
%   * Add date and time to config file
% 
% --------------------------------------------------------------------------
% To Do:
% --------------------------------------------------------------------------
% 
%   * Check restriction values
%   * Figure out the purpose of lines marked with %?%
%   * Check ablation depth in target
%   * Check adsorption coefficient value (Temp and photon energy dependend)
%   * Check crystal binding energy
%   * Check activation energies
%   * Check oxidation energies
%   * Check angular fitting parameter
%   * Check inital maximal velocity
%   * Check velocity delta
%   * Check with data is relevant to print in config file
%   * Improve naming conventions following paper and Matlab guidelines
%   * Correct inital velocity distribution ( f_i(v) )
%   * Implement excitation energy
% 
% --------------------------------------------------------------------------
% Questions:
% --------------------------------------------------------------------------
% 
%   * Why is the crystal binding energy 5x the formation energy?
%   * It says that is 'Oxidation' is set to 1 (true) there is no oxidation
%       possible, is this correct?
%   * What is the target density, and why is it 1?
%   * Why are O2 and Ti2O3 not included in mass array?
%   * Why is only the activation energy of TiO listed, and not of SrO or
%       O2?
%   * What do the three constants in the initial plume settings block do?
%   * Why clear variables from memory within loop, instead of overwriting?
%   * Where does the 0.9 come from in Tom's average velocity?
%   * Why is the average velocity twice divided by the number of atoms in
%       the unit cell?
%   * V_Delta is experimentally fitted according to the paper, but is
%       exactly 5?
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Naming conventions: following MATLAB Style Guidelines v1.3 and loosely
%                       following names used in Wijnand's paper
%   lowers case : variable  (can be altered within script)
%   Camel Case  : Parameter (can be altered by user before running script)
%   UPPER CASE  : CONSTANT  (can not be altered)
%   "t"         : Time
%   "x"         : Distance
%   "rad"       : Angle
%   "n"         : Number of particles
%   "v"         : Velocity
%   "k"         : Number of collisions
%   "p"         : Pressure
%   "temp"      : Temperature
%   "e"         : Energy
%   "mo"        : Moving particles
%   "st"        : Static particles
%   "bg"        : Background gas particles
%   "col"       : Colliding particles
%   "ncol"      : NOT colliding particles 
%   "m"         : Current space matrix of all plume particles
%   "m_1"       : New space matrix of all plume particles
%   "m_bg"      : Current space matrix of all background gas particles
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace
clear
%?%clearvars -except files nn
close all
clc

% Select execution time profiler
%profile off
%?%warning('OFF')

C = PhysicalConstants;

% % Physical constants
% BOLTZMANN = 1.38064852E-23; % Boltzman constant [m2 kg / s2 K]
% AMU = 1.66053906660E-27;    % Atomic mass unit [kg]
% EV  = 1.602176634E-19;      % Electron volt [J]
% PLANCK = 6.62607004E-34;    % Planck constant [m2 kg / s]
% SOL = 299792458;            % Speed of light [m / s]

%% Initial conditions
%--------------------------------------------------------------------------
% File settings
%--------------------------------------------------------------------------

% Save location
directory   = ['C:\Users\Sam\Google Drive\School\Master\Capita selecta' ...
                '\Model_Sam\']; % Parent location
folder      = 'SrTiO3_test';    % Output folder
outputFile  = 'config';         % Output file
%?%SAVE_LOC    = [Dir Folder];  % Save location

% Add comment to output file
outputComment = ('SrTiO3_test');

%--------------------------------------------------------------------------
% Plot settings
%--------------------------------------------------------------------------

% ThreeD      = false;            % Plot 3D results: true / false
% Times_3D    = [5, 15, 20, 30];  % Times for which to generate 3D plots [s]
%?%Arr_sum    = false;            % Calculate arriving particles:  true / false
% Plot_1D     = false;            % Plot 1D results of middle of the plume: true / false
% Plot_BG     = false;            % Show bg in plot: true / false

%--------------------------------------------------------------------------
% Dimensional settings
%--------------------------------------------------------------------------

% Temporal settings
T_0         = 0;        % Start time [s]
T_End       = 18E-6;    % End time [s]
T_Delta     = 1E-7;     % Time step duration [s]

% Spatial settings
X_0         = 0;        % Start position [m]
X_End       = 0.06;     % End position [m]
X_Delta     = 0.5E-4;   % Spatial step size [m]

% Radial settings
Rad_0       = 0;                    % Start angle [deg]
Rad_End     = 90;                   % End angle [deg]
Rad_Delta   = 1;                    % Radial step size [deg]
%?%Rad_N    = 90;                   % Number of total angles
%?%theta = 0 : Rad_Delta : Rad_N;   % Angular array
rad = Rad_0 : Rad_Delta : Rad_End;  % Angular array

% Velocity settings
VEL_MIN   = X_Delta / T_Delta; % Minimal initial velocity
VEL_MAX   = 50000;             % Maximal initial velocity
VEL_DELTA = 5;                 % Velocity step size 0 : (V_MIN/V_Delta) : V_Max;
vel = 0 : (VEL_MIN / VEL_DELTA) : VEL_MAX;
N_VEL = numel(vel);

%--------------------------------------------------------------------------
% Simulation setting
%--------------------------------------------------------------------------

% Restrictions (Check values)
Oxidize             = true;     %?% Enable particle oxidation (true = no oxidization)
Dioxide             = false;    % Let particles form dioxides: true / false
N_Local_Dist_Min    = 1;        % Minimal N per bin
K_Max               = 7;        % Maximum k per particle
N_Min_Rad           = 1E7;      % Minimal N per angle (angles with lower N are not computed)

% Smoothing parameters
Smooth_Fraction = 1E4;  % Fraction to smooth % ???
Smooth_Value    = 15;   % Moving average smooth value (default: 15)

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Atomic masses
% M_A = 87.62  * C.AMU;   % Mass Sr atom [kg]
M_A = 6.94 * C.AMU;     % Mass Li atom [kg]
M_B = 47.867 * C.AMU;   % Mass Ti atom [kg]
M_O = 15.999 * C.AMU;   % Mass O atom [kg]

% Array holding all possible plume compounds
massSpecies = [M_A M_B M_O 2*M_O M_A+M_O M_B+M_O M_B+2*M_O 2*M_B+3*M_O]; % Sr Ti O O2 SrO TiO TiO2 Ti2O3
M_NUMEL = numel(massSpecies);

%?%n_A     = find(M_Array == M_A) - 1;     % place of A atom in the array
%?%n_B     = find(M_Array == M_B) - 1;     % place of B atom in the array
%?%n_O     = find(M_Array == M_O) - 1;     % place of oxygen atom in the array
%?%n_AO    = find(M_Array == M_A+M_O) - 1; % place of AO molecule in the array
%?%n_BO    = find(M_Array == M_B+M_O) - 1; % place of BO molecule in the array
%?%n_BO2   = find(M_Array == M_B+M_O) - 1; % place of BO2 molecule in the array

%?%n_position = [n_A n_B n_O n_AO n_BO n_BO2]; %array of the position in m_p array

% Unit cell
% UC_N_Atoms      = [1 1 3];          % Amount of each atom in STO unit cell [A B O]
nAtomsUC      = [2 1 3];          % Amount of each atom in LTO unit cell [A B O]
UC_N_NUMEL      = numel(nAtomsUC);
UC_Vol          = (3.945E-10)^3;    % Volume of SrTiO3 unit cell [m3]
% UC_E_Formation  = 3.5615 * C.EV;      % Formation energy of STO [J]

% Formation energy of LTO [J]
% https://materialsproject.org/materials/mp-2931/
energyFormation  = 2.990 * C.EV;

% Optical properties
adsorptionC = 0.77; % Adsorbtion coefficient SrTiO3 (Temp and laser energy dependend) (check value)

% Target
%?%Density_Target = 1;  % Target density

% Energies
E_TiO = 0.117 * C.EV;                 % Activation energy for forming TiO [J]
E_O_low = 0.5 * C.EV;                 % Lower bound oxidation energy [J]
E_O_high_Sr = 4.418 * C.EV; %* Oxidize; % Upper bound oxidation energy Sr [J]
E_O_high_Ti = 2.276 * C.EV; %* Oxidize; % Upper bound oxidation energy Ti [J]

heatTarget = 0;  % No significant heat loss into the target (ceramic)
energyExcitation = zeros(1, M_NUMEL);   % No significant excitation of species

%?%t_p = 4;     %number of start particles +1
%?%% t_o = 4;   %number of oxidation states
%?%% t_oo = 2;  %placement oxide particle.
%?%% t_pp = 8;  %total number of particles

% Atomic radii [pm]
% Data from: doi.org/10.1063/1.1712084
% Assume a molecule has atomic radii of the sum of its components
Sr_R = 219;
Ti_R = 176;
O_R = 48;

% Collision cross section array of plume species with background gas molecule [m2]
Sigma_Bg = (([Sr_R ...              % Sr-O2
            Ti_R ...                % Ti-O2
            O_R ...                 % O-O2
            (Sr_R + O_R) ...        % SrO-O2
            (Ti_R + O_R) ...        % TiO-O2
            (Ti_R + 2*O_R) ...      % TiO2-O2
            ] + 2*O_R).*(10^-12).^2) .* pi;

%?% Opp_plasma = ([(20E-11 + 4.8E-11*2) ...
%    (17.6E-11 + 4.8E-11*2) ...
%    (3*4.8E-11) ...
%    (20E-11 + 4.8E-11*3) ...
%    (17.6E-11 + 4.8E-11*3) ...
%    (17.6E-11 + 4.8E-11*4)] .^2) .* pi; %?% [Sr-O2, Ti-O2, O-O2, SrO-O2, TiO-O2, TiO2-O2] [m^2]

%--------------------------------------------------------------------------
% Deposition settings
%--------------------------------------------------------------------------

% Background gas
P_Bg    = 0.1;  % Pressure of background gas during depostion [mbar]
Temp_Bg = 300;  % Temperature of background gas [K]

% System dimensions
X_Ts = 0.05;    % Target-substrate distance [m]

% Laser parameters
energyLaser = 0.0299;   % Laser energy [J]
Spot_Width  = 1E-3;     % Laser spot width [m]
Spot_Height = 2E-3;     % Laser spot height [m]
Spot_Depth  = 100E-9;   % Ablation depth in the target [m]

%--------------------------------------------------------------------------
% Initial plume settings
%--------------------------------------------------------------------------

Rad_Fit = 16;           % Radial sharpe of the plasma
velDisWidth = 1654;         %?% Width velocity distribution
%?%dis_NPoints =  100;  % Number of datapoints in the initial velocity distribution

%--------------------------------------------------------------------------
%% Constant calculation based on parameters

% Background gas
P_BG_PASCAL = 100 * P_Bg;                       % Convert pressure in Pa (1 mbar = 10^2 Pa)
DENSITY_BG  = P_BG_PASCAL / (C.BOLTZMANN * Temp_Bg);    % Density N/V of background gas (ideal gas law)

% Laser
SPOT_VOLUME = Spot_Width * Spot_Height * Spot_Depth;    % Laser spot volume [m3]

% Unit cell
UC_N_TOTAL      = sum(nAtomsUC);              % Total number of atoms in unit cell
UC_E_BINDING    = UC_N_TOTAL * energyFormation;  % Binding energy crystal (5 times formation energy per atom) [eV]

%--------------------------------------------------------------------------
%% Initial ditribution
nUCAblated = SPOT_VOLUME / UC_Vol;  % Number of ablated unit cells
% N_UC = (SPOT_VOLUME * UC_N_TOTAL) / UC_Vol; % * Density_Target)); % Total number of ablated atoms per UC
%?%N_TOTAL = (UC_N_Atoms * N_UC) / UC_N_TOTAL;  % Total number of ablated atoms

n_atom_temp = zeros(numel(rad) - 1, 1); % Pre-allocate memory

% Compute initial partical distribution
for irad = 1 : (numel(rad) - 1)
    n_atom_temp(irad) = ((4/3)*pi * X_Delta^3) .* (cosd(rad(irad)) ...
        - cosd(rad(irad + 1))) .* cosd(rad(irad)).^Rad_Fit;
end

% Number of atoms per angle
n_atom_rad = (n_atom_temp .* nUCAblated) ./ sum(n_atom_temp); %?%

%--------------------------------------------------------------------------
%% Write settings into config file

% Make directory if none exists
if ~isfolder([directory folder])
    mkdir([directory folder])
end

% Generate file name
file_name = [directory folder '\' outputFile '_' datestr(now, 'yyyy-mm-dd_HHMM') '.txt'];

% Open file for writing
fid = fopen(file_name, 'wt');

fprintf(fid, ['\tCONFIGURATION FILE\t\t' datestr(now, 'yyyy-mm-dd\tHH:MM') '\n\n' ...
    'Temporal settings [s]:' '\t\t\t' num2str(T_0) ' : ' num2str(T_Delta) ' : ' num2str(T_End) '\n' ...
    'Spatial settings [m]:' '\t\t\t' num2str(X_0) ' : ' num2str(X_Delta) ' : ' num2str(X_End) '\n' ...
    'Angular settings [deg]:' '\t\t\t' num2str(Rad_0) ' : ' num2str(Rad_Delta) ' : ' num2str(Rad_End) '\n'...
    'Velocity limits [m/s]:' '\t\t\t' 'min ' num2str(VEL_MIN) ', max ' num2str(VEL_MAX) ', dv ' num2str(VEL_DELTA) '\n' ...
    'Background pressue [mbar]:' '\t\t' num2str(P_Bg) '\n' ...
    'Oxidation [true/false]:' '\t\t\t' mat2str(Oxidize) '\n' ...
    'Collision cross-sections [m2]:' '\t\t' mat2str(Sigma_Bg,2) '\t' '[Sr-O2, Ti-O2, O-O2, SrO-O2, TiO-O2, TiO2-O2]' '\n' ...
    '\n\tComments:' '\n' outputComment]);

%{
fprintf(fid, ['Number of type of particles,\t\t' num2str(t_p)  '\n']);
fprintf(fid, ['UC volume,\t\t' num2str(UC_vol)  '\n']);
fprintf(fid, ['Mass atom,\t\t' num2str(m_B)  '\n']);
fprintf(fid, ['Mass Bg,\t\t' num2str(m_O)  '\n']);
fprintf(fid, ['Mass particles,\t\t' num2str(m_P)  '\n']);
fprintf(fid, ['Laser Energy,\t\t' num2str(dis_E_x)  '\n']);
fprintf(fid, ['N_points,\t' num2str(dis_NPoints)  '\n']);
fprintf(fid, ['absorption coefficient,\t\t' num2str(A_coeff)  '\n']);
fprintf(fid, ['binding energy crystal,\t\t' num2str(Eb_UC)  '\n']);
fprintf(fid, ['Density target,\t\t' num2str(Density_target)  '\n']);
fprintf(fid, ['Number of atoms,\t' num2str(N_atoom_totaal)  '\n']);
fprintf(fid, ['CosN,\t\t\t' num2str(cos_theta)  '\n']);
fprintf(fid, ['Degs,\t\t\t' num2str(N_rad)  '\n']);
fprintf(fid, ['Velocity distribution,\t' num2str(sqrt((dis_E_x/(A_coeff*N_atoom_totaal/sum(T_sum))-Eb_UC)*(2/(sum(T_sum)*m_P(n_A+1)))))  '\n']);
fprintf(fid, ['Velocity Width,\t\t' num2str(V_delta)  '\n']);
fprintf(fid, ['E_O_upper,\t' num2str(E_O_high_Sr)  '\n']);
fprintf(fid, ['E_O_lower,\t\t' num2str(E_O_low)  '\n']);
fprintf(fid, ['Spot size,\t\t' num2str(Spot_size)  '\n']);
fprintf(fid, ['3D plot dried?,\t\t' num2str(dried)  '\n']);
fprintf(fid, ['MO2 possible,\t\t' num2str(MO2)  '\n']);
fprintf(fid, ['Vminimaal,\t\t' num2str(V_minimaal)  '\n']);
fprintf(fid, ['Volume BG,\t\t' num2str(4 * pi/3 * (((1)*X_delta)^3-((1-1)*X_delta)^3) * (-cosd(Rad_delta)+cosd(Rad_delta-Rad_delta)))  '\n']);
fprintf(fid, ['Target substrate distance,\t\t' num2str(X_TS)  '\n']);
fprintf(fid, ['maximum number of collisions,\t\t' num2str(limbots)  '\n']);
fprintf(fid, ['Minimal amount of particles,\t\t' num2str(N_local_dist_min)  '\n']);
fprintf(fid, ['smooth_fraction,\t\t' num2str(smooth_fraction)  '\n']);
fprintf(fid, ['smooth_value,\t\t' num2str(smooth_value_set)  '\n']);
%}

fclose(fid);

%--------------------------------------------------------------------------
%% Pre-allocations

%--------------------------------------------------------------------------
%% Calculate the initial velocity praticle distribution
nVelDisInit = initialVelocityDistribution( vel, velDisWidth, ...
    nUCAblated, nAtomsUC, massSpecies, energyLaser, energyFormation, ...
    energyExcitation, heatTarget, adsorptionC );

%--------------------------------------------------------------------------
%% Main program

for theta = 1 : 1 %(numel(rad) - 1)
    % Initiate
    n_atom = 1;
    
    % Number of atoms at current angle
    if n_atom_rad(theta) > N_Min_Rad
        n_atom = n_atom_rad(theta);
    end
end
