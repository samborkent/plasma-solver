%% PLASMA SOLVER v5.0.0 23-02-21
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
%   dist    : Distribution
%   init    : Initial
%   n       : Number of particles
%   bg      : Background gas
%   uc      : Unit cell
%   bin     : Computational bin with similar position, velocity, etc.
%
% Field names and corresponding index value (see Field enumaration class):
%   Field.veloBins    (1)
%        .nParticles  (2)
%        .nCollisions (3)
%        .species     (4)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace
clc, clear, close all

%% File settings

% Get current path
currentPath = pwd;

% Add 
addpath([currentPath '/functions/']);
addpath([currentPath '/classes/']);

% Save location
directory   = currentPath;      % Parent location
folder      = 'results';        % Output folder
fileName    = 'config';         % Output file

% Add comment to output file
commentString = [...
    'Type whatever you want here. For example, this script was the ' ...
    'worst sturctured code I have ever seen.' ...
    ];

%% Plot settings
createConfigBool    = true;    % Create a configuration file
plotVelDisInit      = false;    % Plot initial velocity distribution
saveVelDisInit      = false;    % Save the initial velocity distribution plot

%% Create class instances holding constants and material properties
C   = PhysicalConstants;
PT  = PeriodicTable;
UC  = UnitCells;

%% Initial conditions
%--------------------------------------------------------------------------
% Limits
%--------------------------------------------------------------------------

nParticleMin = 1;    % Minimal number of particles per bin

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 18E-6;    % End time [s]
% timeDelta   = 1E-7;     % Time step duration [s]
timeDelta   = 1E-6;

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
% radiusDelta = 0.5E-4;   % Spatial step size [m]
radiusDelta = 0.001;

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Velocity limits
veloMin         = 0;                        % Minimal initial velocity
veloMax         = 5E4;                      % Maximal initial velocity
veloDelta       = radiusDelta / timeDelta;  % Velocity delta to move to next velocity bin
veloStepsize    = 5;                        % Velocity step size (veloMin : (veloDelta/veloStepsize) : V_Max)

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Unit cell
uc              = UC.Li4Ti5O12;  
atomUC          = uc.ELEMENTS;      % Atoms in unit cell
nAtomUC         = uc.AMOUNT;        % Amount of each atom in unit cell [A B O]
nAtomUCNumel    = numel(nAtomUC);
nAtomUCSum      = sum(nAtomUC);
ucVolume        = uc.VOLUME;        % Volume of unit cell [m^3]

% Select A and B atoms  
massA = atomUC(1).MASS;     % Mass A atom [kg]
massB = atomUC(2).MASS;     % Mass B atom [kg]
massO = atomUC(end).MASS;   % Mass O atom [kg]

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
radiusO     = atomUC(end).MASS; % Oxygen atom
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
absorption          = 0.66; % Adsorption coefficient SrTiO3 (check value)
heatTarget          = 0;    % No significant heat loss into the target (ceramic)

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

% Background gas density (ideal gas law) 
bgDensity = bgPressure / (C.BOLTZMANN * bgTemperature);

% Ablation volume [m^3]
ablationVolume  = spotWidth * spotHeight * ablationDepth;

% Laser energy [J]
% energyLaser = 0.0520; % Tom
energyLaser = (laserFluence * 10^4) * (spotWidth * spotHeight);

% Unit cell
energyBinding = energyFormation;    % Binding energy crystal [J]

% Dimensional axis
time    = timeMin   : timeDelta                 : timeMax;   % Temporal axis
radius  = radiusMin : radiusDelta               : radiusMax; % Radial axis
angle   = angleMin  : angleDelta                : angleMax;  % Angular axis
velo    = veloMin   : veloDelta / veloStepsize  : veloMax;   % Velocity array
veloBins = veloMin : veloDelta : veloMax - veloDelta;

% Dimensional sizes
nTime   = numel(time);
nRadius = numel(radius);
nAngle  = numel(angle);
nVelo   = numel(velo);
nVeloBins = numel(veloBins);

%--------------------------------------------------------------------------
%% Initial ditribution
% Number of ablated unit cells
% nUCAblated = (ablationVolume / ucVolume) / nAtomUCSum; % Tom
nUCAblated      = ablationVolume / ucVolume;    % Number of ablated unit cells
nAtomAblated    = nUCAblated * nAtomUCSum;      % Number of ablated atoms

% Pre-allocate memory
nParticleDistribution = zeros(1, nAngle - 1);

% Compute initial angular particle distribution (eq. 5)
for iAngle = 1 : (nAngle - 1)
    nParticleDistribution(iAngle) = ((4/3)*pi * radiusDelta^3) ...
        .* abs( cosd(angle(iAngle + 1)) - cosd(angle(iAngle)) ) ...
        .* cosd(angle(iAngle)).^cosPowerFit; %?%
end

% Normalize and multiply by the total number of ablated atoms
nParticleDistribution = (nParticleDistribution .* nAtomAblated) ./ sum(nParticleDistribution);

%--------------------------------------------------------------------------
%% Write settings into config file
if createConfigBool
    createConfigFile( directory, folder, fileName, commentString, time, ...
        radius, angle, velo, bgPressure );
end

%--------------------------------------------------------------------------
%% Preallocation
% Number of property fields per velocity bin
nFieldsPerVeloBin = numel(enumeration('Field'));

% Plasma plume particle matrix
plasmaMatrix = zeros(nVeloBins, nFieldsPerVeloBin, nRadius - 1);

% Insert velocity bins
for iRadius = 1 : (nRadius - 1)
    plasmaMatrix(:, Field.veloBins, iRadius) = veloBins;
end

% Copy plasma matrix
bgMatrix = plasmaMatrix;

for iAngle = 1 : 1
    %% Main program per angle
    % Initiate number of atoms
    nParticleAngle = nParticleMin;
    
    % Number of atoms at current angle
    if nParticleDistribution(iAngle) > nParticleMin
        nParticleAngle = nParticleDistribution(iAngle);
    end
    
    % Calculate the initial particle velocity distribution
    veloDistInit = initialVelocityDistribution( plotVelDisInit, ...
        saveVelDisInit, veloBins, veloDistributionWidth, nUCAblated, uc, ...
        energyLaser, energyBinding, heatTarget, absorption, nParticleAngle );
    
    % Only plot and save initial particle velocity distribution for the
    % center of the plume
    if iAngle == 1
        plotVelDisInit = false;
        saveVelDitInit = false;
    end
    
    plasmaMatrix(:, Field.nParticles, 1) = veloDistInit;
    
    for iTime = 1 : nTime
        %% Main program per timestep
        
        for iRadius = 1 : (nRadius - 1)
            %% Main program per radial bin
            
            % Calculate bin volumes and background particle density per bin
            binVolume = (4/3)*pi ...
                * ( (radius(iRadius + 1))^3 ...
                - (radius(iRadius))^3 ) ...
                * ( cosd(angle(iAngle)) ...
                - cosd(angle(iAngle + 1)) ); % Bin volume [m^3]

            % Background particle density per bin
            bgBinNParticle = bgDensity * binVolume;

            % Insert in background particle matrix
            bgMatrix(:, Field.nParticles, iRadius) = bgBinNParticle;
            
            % Loop though velocity bins from highest to lowest
            for iVelo = numel(veloBins) : -1 : 1
                %% Main program per velocity bin
                
                % Store matrix properties in variables
                veloPlasma  = plasmaMatrix(iVelo, Field.veloBins, iRadius);
                nPlasma     = plasmaMatrix(iVelo, Field.nParticles, iRadius);
                veloBg      = bgMatrix(iVelo, Field.veloBins, iRadius);
                nBg         = bgMatrix(iVelo, Field.nParticles, iRadius);
                
                % Average distance traveled in this time step
                radiusTraveled = abs(veloPlasma * timeDelta);
                
                % Number of radius bins traveles
                nRadiusBinTraveled = radiusTraveled / radiusDelta;
                
                % Limit number of traveled radius bins to prevent index
                % out-of-bounds error
                if nRadiusBinTraveled > (nRadius - 2)
                    nRadius = nRadius - 2;
                end
                
                % If the bin has not enough particles
                if (nPlasma < nParticleMin) && (nRadiusBinTraveled < 1)
                    colChanceStatic     = zeros(1, nRadiusBinTraveled);
                    colChanceKinetic    = zeros(1, nRadiusBinTraveled);
                % If bin has enough particles
                else
                    % Loop through distance traveled in steps radiusDelta
                    for iRadiusDelta = iRadius : nRadiusBinTraveled
                        % Calculate collision probability per radius bin
                        
                    end % For radiusDelta
                end % If number of particles
            end % For velocity
        end % For radius
    end % For time
end % For angle
