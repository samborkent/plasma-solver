%% PLASMA SOLVER v5.0.0 17-02-21
%
% Authors:   Tom Wijnands, Bart Boonsma, Sam Borkent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Naming conventions:
% Loosly resemble paper names and follow MATLAB style guidelines v1.3
%
%   radius  : Radius / Distance
%   velo    : Velocity
%   n       : Number of particles
%   bg      : Background gas
%   uc      : Unit cell
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace
clc, clear, close all

CONSTANT = PhysicalConstants;

%% Initial conditions
%--------------------------------------------------------------------------
% File settings
%--------------------------------------------------------------------------

% Save location
directory   = ['C:\Users\Sam\Google Drive\School\Master\Capita selecta' ...
                '\Model_Sam\']; % Parent location
folder      = 'SrTiO3_test';    % Output folder
fileName    = 'config';         % Output file

% Add comment to output file
commentString = ('SrTiO3_test');

%--------------------------------------------------------------------------
% Plot settings
%--------------------------------------------------------------------------

plotVelDisInit  = true;           % Plot initial velocity distribution (true / false)

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 18E-6;    % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 0.5E-4;   % Spatial step size [m]

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 1;        % Radial step size [deg]

% Velocity limits
veloMin     = radiusDelta / timeDelta;  % Minimal initial velocity
veloMax     = 6E4;                      % Maximal initial velocity
veloDelta   = 5;                        % Velocity step size 0 : (V_MIN/V_Delta) : V_Max;

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Select A and B atoms
massA = CONSTANT.MASS_Li;   % Mass Li atom [kg]
massB = CONSTANT.MASS_Ti;   % Mass Ti atom [kg]
massO = CONSTANT.MASS_O;    % Mass O atom [kg]

% Array holding all possible plume compounds
massSpecies = [ massA ...               % atomic lithium        [Li +1]
                massB ...               % atomic titanium       [Ti +4]
                massO ...               % atomic oxygen         [O -2]
                2*massA ...             % dilithium             [Li2 +2]
                2*massO ...             % molecular oxygen      [O2]
                massA + massO ...       % lithium monoxide      [LiO -1]
                2*massA + massO ...     % lithium oxide         [Li2O]
                2*massA + 2*massO ...   % lithium peroxide      [Li2O2 -2]
                4*massO ...             % tetraoxygen           [O4 -4]
                massB + massO ...       % titanium(II) oxide    [TiO +2]
                massB + 2*massO ...     % titanium dioxide      [TiO2]
                2*massB + massO ...     % dititanium monooxide  [Ti2O +6]
                2*massB + 3*massO ...   % titanium(III) oxide   [Ti2O3 +2]
                3*massB + massO         % trititanium oxide     [Ti3O +10]
                ];
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
energyFormation  = 2.990 * CONSTANT.EV;

% Optical properties
adsorptionC = 0.77; % Adsorbtion coefficient SrTiO3 (Temp and laser energy dependend) (check value)

% Target
%?%Density_Target = 1;  % Target density

% Energies
E_TiO = 0.117 * CONSTANT.EV;                 % Activation energy for forming TiO [J]
E_O_low = 0.5 * CONSTANT.EV;                 % Lower bound oxidation energy [J]
E_O_high_Sr = 4.418 * CONSTANT.EV; %* Oxidize; % Upper bound oxidation energy Sr [J]
E_O_high_Ti = 2.276 * CONSTANT.EV; %* Oxidize; % Upper bound oxidation energy Ti [J]

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
pressureBG = 0.1;  % Pressure of background gas during depostion [mbar]
temperatureBG = 300;  % Temperature of background gas [K]

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
P_BG_PASCAL = 100 * pressureBG;                       % Convert pressure in Pa (1 mbar = 10^2 Pa)
DENSITY_BG  = P_BG_PASCAL / (CONSTANT.BOLTZMANN * temperatureBG);    % Density N/V of background gas (ideal gas law)

% Laser
SPOT_VOLUME = Spot_Width * Spot_Height * Spot_Depth;    % Laser spot volume [m3]

% Unit cell
UC_N_TOTAL      = sum(nAtomsUC);              % Total number of atoms in unit cell
UC_E_BINDING    = UC_N_TOTAL * energyFormation;  % Binding energy crystal (5 times formation energy per atom) [eV]

time = timeMin : timeDelta : timeMax;           % Temporal axis
radius = radiusMin : radiusDelta : radiusMax;   % Radial axis
angle = angleMin : angleDelta : angleMax;       % Angular axis
velo = 0 : (veloMin / veloDelta) : veloMax;     % Velocity array

nTime = numel(time);
nRadius = numel(radius);
nAngle = numel(angle);
nVelo = numel(velo);

%--------------------------------------------------------------------------
%% Initial ditribution
nUCAblated = SPOT_VOLUME / UC_Vol;  % Number of ablated unit cells
% N_UC = (SPOT_VOLUME * UC_N_TOTAL) / UC_Vol; % * Density_Target)); % Total number of ablated atoms per UC
%?%N_TOTAL = (UC_N_Atoms * N_UC) / UC_N_TOTAL;  % Total number of ablated atoms

n_atom_temp = zeros(nAngle - 1, 1); % Pre-allocate memory

% Compute initial partical distribution
for iAngle = 1 : (nAngle - 1)
    n_atom_temp(iAngle) = ((4/3)*pi * radiusDelta^3) .* (cosd(angle(iAngle)) ...
        - cosd(angle(iAngle + 1))) .* cosd(angle(iAngle)).^Rad_Fit;
end

% Number of atoms per angle
n_atom_rad = (n_atom_temp .* nUCAblated) ./ sum(n_atom_temp); %?%

%--------------------------------------------------------------------------
%% Write settings into config file

createConfigFile( directory, folder, fileName, commentString, time, ...
    radius, angle, velo, pressureBG );

%--------------------------------------------------------------------------
%% Calculate the initial particle velocity distribution
nVelDisInit = initialVelocityDistribution( plotVelDisInit, velo, ...
    velDisWidth, nUCAblated, nAtomsUC, massSpecies, energyLaser, ...
    energyFormation, energyExcitation, heatTarget, adsorptionC );

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
