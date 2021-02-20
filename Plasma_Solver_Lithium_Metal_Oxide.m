%% PLASMA SOLVER v5.0.0 17-02-21
%
% Authors:   Tom Wijnands, Bart Boonsma, Sam Borkent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Naming conventions:
%   Loosly resemble paper names and follow MATLAB style guidelines v1.3
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
PT = PeriodicTable;
UC = UnitCells;

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
commentString = [...
    'Type whatever you want here. For example, this code was one ' ...
    'of the worst sturctured code I have ever seen.' ...
    ];

%--------------------------------------------------------------------------
% Plot settings
%--------------------------------------------------------------------------

plotVelDisInit  = true; % Plot initial velocity distribution (true / false)

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
angleDelta  = 3;        % Radial step size [deg]

% Velocity limits
veloMin     = radiusDelta / timeDelta;  % Minimal initial velocity
veloMax     = 20E4;                      % Maximal initial velocity
veloDelta   = 100;
% veloDelta   = 5;                        % Velocity step size 0 : (V_MIN/V_Delta) : V_Max;

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Unit cell
uc              = UC.Li4Ti5O12;
atomUC          = [PT.Li PT.Ti PT.O];   % Atoms in unit cell
nAtomUC         = [4 5 12];              % Amount of each atom in unit cell [A B O]
nAtomUCNumel    = numel(nAtomUC);
nAtomUCSum      = sum(nAtomUC);
ucVolume        = uc.VOLUME;            % Volume of unit cell [m^3]

% Select A and B atoms  
massA = atomUC(1).MASS;     % Mass A atom [kg]
massB = atomUC(2).MASS;     % Mass B atom [kg]
massO = atomUC(3).MASS;     % Mass O atom [kg]

% Array holding all possible plume compounds
massArray = [ massA ...               % atomic lithium        [Li +1]
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
              3*massB + massO  ...    % trititanium oxide     [Ti3O +10]
              3*massB + 5*massO       % trititanium pentoxide [Ti3O5 +2]
            ];
nSpecies = numel(massArray);

% Atomic radii [m]
% Assume a molecule has atomic radii of the sum of its components
radiusA     = atomUC(1).MASS;   % A atom: lithium
radiusB     = atomUC(2).MASS;   % B atom: titanium
radiusO     = atomUC(3).MASS;   % Oxygen atom
radiusBG    = 2*PT.O.RADIUS;    % Background gas molecule

% Collision cross section array of plume species with background gas molecule [m2]
Sigma_Bg = (([ radiusA ...                  % A-O2
               radiusB ...                  % B-O2
               radiusO ...                  % O-O2
               (radiusA + radiusO) ...      % AO-O2
               (radiusB + radiusO) ...      % BO-O2
               (radiusB + 2*radiusO) ...    % BO2-O2
             ] + radiusBG).^2) .* pi;

% Formation energy of monoclinic Li4Ti5O12 [J] (get spinel value)
energyFormation  = uc.ENERGY_FORMATION;

% Target properties
adsorptionC         = 0.66;                    % Adsorption coefficient SrTiO3 (check value)
heatTarget          = 0;                       % No significant heat loss into the target (ceramic)

%--------------------------------------------------------------------------
% Deposition settings
%--------------------------------------------------------------------------

% Background gas
bgPressure      = 0.1E2;    % Background gaspressure during depostion [Pa]
bgTemperature   = 300;      % Temperature of background gas [K]

% Laser parameters
spotWidth       = 1E-3;     % Laser spot width [m]
spotHeight      = 2.1E-3;   % Laser spot height [m]
ablationDepth   = 100E-9;   % Ablation depth in the target [m]
laserFluence    = 2.0;      % Laser fluence [J / cm^2]   

%--------------------------------------------------------------------------
% Initial plume settings
%--------------------------------------------------------------------------

cosPowerFit = 25;               % Radial sharpe of the plasma (check value)
veloDistributionWidth = 1500;   % Initial particle velocity distribution width (check value)

%--------------------------------------------------------------------------
%% Constant calculation based on parameters

% Background gas density (ideal gas law) []
bgDensity = bgPressure / (CONSTANT.BOLTZMANN * bgTemperature);

% Ablation volume [m3]
ablationVolume  = spotWidth * spotHeight * ablationDepth;

% Laser energy [J]
% energyLaser = 0.0520; % Tom
energyLaser = (laserFluence * 10^4) * (spotWidth * spotHeight);

% Unit cell
energyBinding = energyFormation;    % Binding energy crystal [J]

% Dimensional axis
time    = timeMin               : timeDelta     : timeMax;      % Temporal axis
radius  = radiusMin             : radiusDelta   : radiusMax;    % Radial axis
angle   = angleMin + angleDelta : angleDelta    : angleMax;     % Angular axis
velo    = veloMin               : veloDelta     : veloMax;      % Velocity array
% velo  = 0 : (veloMin / veloDelta) : veloMax;

% Dimensional sizes
nTime   = numel(time);
nRadius = numel(radius);
nAngle  = numel(angle);
nVelo   = numel(velo);

%--------------------------------------------------------------------------
%% Initial ditribution
% Number of ablated unit cells
% nUCAblated = (ablationVolume / ucVolume) / nAtomUCSum; % Tom
nUCAblated = ablationVolume / ucVolume; % Number of ablated unit cells
nAtomAblated = nUCAblated * nAtomUCSum; % Number of ablated atoms

% Pre-allocate memory
nParticleAngle = zeros(1, nAngle - 1);

% Compute initial angular particle distribution
for iAngle = 1 : (nAngle - 1)
    nParticleAngle(iAngle) = ((4/3)*pi * radiusDelta^3) .* ...
        (cosd(angle(iAngle)) - cosd(angle(iAngle + 1))) .* ...
        cosd(angle(iAngle)).^cosPowerFit; %?%
end

% Normalize and multiply by the total number of ablated atoms
nParticleAngle = (nParticleAngle .* nAtomAblated) ...
                    ./ sum(nParticleAngle);

%--------------------------------------------------------------------------
%% Write settings into config file
createConfigFile( directory, folder, fileName, commentString, time, ...
    radius, angle, velo, bgPressure );

%--------------------------------------------------------------------------
%% Calculate the initial particle velocity distribution
[nVeloDistributionInitial, E_k] = initialVelocityDistribution( plotVelDisInit, ...
    velo, veloDistributionWidth, nUCAblated, atomUC, nAtomUC, ...
    energyLaser, energyBinding, heatTarget, adsorptionC, nParticleAngle(2) );

%--------------------------------------------------------------------------
%% Main program

E_k = E_k / CONSTANT.EV;

% for theta = 1 : 1 %(numel(rad) - 1)
%     % Initiate
%     n_atom = 1;
%     
%     % Number of atoms at current angle
%     if n_atom_rad(theta) > N_Min_Rad
%         n_atom = n_atom_rad(theta);
%     end
% end
