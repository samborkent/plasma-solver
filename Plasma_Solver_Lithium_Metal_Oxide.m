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
%
%   Field names and corresponding index value (see Field enumaration class):
%       Field.nParticles    (1)
%            .nCollisions   (2)
%
%   Example:
%       nParticleVelo: number of plasma plume particles per velocity bin
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
% directory   = currentPath;    % Parent location
folder      = 'results';        % Output folder
fileName    = 'config';         % File name of configuration file

% Add comment to output file
commentString = [...
    'Type whatever you want here. For example, this script was the ' ...
    'worst structured code I have seen.' ...
    ];

%% Plot settings
createConfigBool    = false;    % Create a configuration file
plotVelDisInit      = true;    % Plot initial velocity distribution
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
timeDelta   = 1E-7;     % Time step duration [s]
% timeDelta   = 1E-6;

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
% radiusDelta = 0.5E-4;   % Spatial step size [m]
radiusDelta = 0.0005;

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Velocity limits
veloMin         = 0;                        % Minimal initial velocity
veloMax         = 4E4;                      % Maximal initial velocity
veloDelta       = radiusDelta / timeDelta;  % Velocity delta to move to next velocity bin
veloStepsize    = 5;                        % Velocity step size (veloMin : (veloDelta/veloStepsize) : V_Max)

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Specified target density [kg / m^3]
% Set to 0 for single crystal target
% targetDensity   = 2910;
targetDensity = 0;

% Unit cell ( start with UC. )
% Options: TiO2, SrTiO3, Li2O, LiO2, LiTi2O4, Li2TiO3, Li4Ti5O12
% uc = [UC.Li4Ti5O12 UC.SrTiO3];    % Array holding unit cells in target
uc = UC.TiO2;

% Mixture ratio of species in target
ucRatio           = [1 3];

% Set ratio to 1 for one component target
if numel(uc) == 1
	ucRatio = 1;
end

% Check if number of given materials in target if the same as the number of
% elements in the mixture ratio array
if numel(uc) ~= numel(ucRatio)
    error(['Number of elements of mixture ratio array must be the ' ...
           'same as the unit cell array.']);
end

ucNumel         = numel(ucRatio);

% % Initialize atom arrays
atomUC          = [];
atomUCTemp      = [];
nAtomUCTemp     = [];
ucVolume        = zeros(1, ucNumel);

% Fill atom arrays
for i = 1 : numel(uc)
    atomUCTemp  = [atomUCTemp   uc(i).ELEMENTS];
    nAtomUCTemp = [nAtomUCTemp  (uc(i).AMOUNT .* ucRatio(i))];
    ucVolume(i) = uc(i).VOLUME;
end

% Normalize atom amount array
nAtomUCTemp = nAtomUCTemp ./ sum(ucRatio);
nAtomUCTemp = nAtomUCTemp ./ min(nAtomUCTemp);

% Initialize atom number array
elementNumberTemp = zeros(1, numel(atomUCTemp));

% Fill atom number array
for i = 1 : numel(atomUCTemp)
    elementNumberTemp(i) = atomUCTemp(i).NUMBER;
end

% Get sorted unique atom array
[elementNumber, elementIndex, ~] = unique(elementNumberTemp, 'sorted');

% Check if there is oxygen
if ~isempty(elementNumber(elementNumber == 8))
    % Place oxygen at the end of the atom array
    elementIndex(end + 1)               = elementIndex(elementNumber == 8);
    % Remove from old array
    elementIndex(elementNumber == 8)    = [];
end

% Initialize amount of atom array
nAtomUC = zeros(1, numel(elementNumber));

% Look through unique atoms
for i = 1 : numel(elementNumber)
    % Insert unique atoms in array
    atomUC      = [atomUC atomUCTemp(elementIndex(i))];
    % Insert the calculated ratio of each atom in atom amount array
    nAtomUC(i)  = sum(nAtomUCTemp(elementNumberTemp == atomUC(i).NUMBER));
end

nAtomUCNumel    = numel(nAtomUC);
nAtomUCSum      = sum(nAtomUC);

% Single crystal density as weighted average of the components in the target [kg / m^3]
scDensity = 0;
for i = 1 : numel(uc)
    scDensity = scDensity + (uc(i).DENSITY * ucRatio(i));
end
scDensity = scDensity / sum(ucRatio);

% Calculate target density ratio compared to single crystal density
if targetDensity
    densityRatio = targetDensity / scDensity;
else
    densityRatio = 1;
end

% Mass of background particles
bgMass = 2*PT.O.MASS;

% Select A and B atoms  
massA = atomUC(1).MASS;     % Mass A atom [kg]
massB = atomUC(2).MASS;     % Mass B atom [kg]
massO = atomUC(end).MASS;   % Mass O atom [kg]

% % Array holding all possible plume compounds
% massArray = [ massA ...               % atomic lithium        [Li +1]
%               massB ...               % atomic titanium       [Ti +4]
%               massO ...               % atomic oxygen         [O -2]
%               2*massA ...             % dilithium             [Li2 +2]
%               2*massO ...             % molecular oxygen      [O2]
%               massA + massO ...       % lithium monoxide      [LiO -1]
%               2*massA + massO ...     % lithium oxide         [Li2O]
%               2*massA + 2*massO ...   % lithium peroxide      [Li2O2 -2]
%               4*massO ...             % tetraoxygen           [O4 -4]
%               massB + massO ...       % titanium(II) oxide    [TiO +2]
%               massB + 2*massO ...     % titanium dioxide      [TiO2]
%               2*massB + massO ...     % dititanium monooxide  [Ti2O +6]
%               2*massB + 3*massO ...   % titanium(III) oxide   [Ti2O3 +2]
%               3*massB + massO  ...    % trititanium oxide     [Ti3O +10]
%               3*massB + 5*massO       % trititanium pentoxide [Ti3O5 +2]
%             ];
% nSpecies = numel(massArray);
% 
% % Atomic radii [m]
% % Assume a molecule has atomic radii of the sum of its components
% radiusA     = atomUC(1).MASS;   % A atom: lithium
% radiusB     = atomUC(2).MASS;   % B atom: titanium
% radiusO     = atomUC(end).MASS; % Oxygen atom
% radiusBG    = 2*PT.O.RADIUS;    % Background gas molecule
% 
% % Collision cross section array of plume species with background gas molecule [m2]
% Sigma_Bg = (([ radiusA ...                  % A-O2
%                radiusB ...                  % B-O2
%                radiusO ...                  % O-O2
%                (radiusA + radiusO) ...      % AO-O2
%                (radiusB + radiusO) ...      % BO-O2
%                (radiusB + 2*radiusO) ...    % BO2-O2
%              ] + radiusBG).^2) .* pi;

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
% energyBinding = energyFormation;    % Binding energy crystal [J]

% Dimensional axis
angle   = angleMin  : angleDelta                : angleMax - angleDelta;   % Angular axis
time    = timeMin   : timeDelta                 : timeMax - timeDelta;     % Temporal axis
radius  = radiusMin : radiusDelta               : radiusMax;               % Radial axis
velo    = veloMin   : veloDelta / veloStepsize  : veloMax - veloDelta;     % Velocity array

% Dimensional sizes
nTime   = numel(time);
nRadius = numel(radius);
nAngle  = numel(angle);
nVelo   = numel(velo);

% Number of radius bins traveled on average per time step per velocity bin
nRadiusBinTraveled = round((velo .* timeDelta) ./ radiusDelta);

%--------------------------------------------------------------------------
%% Initial ditribution
% Number of ablated unit cells
% nUCAblated = (ablationVolume / ucVolume) / nAtomUCSum; % Tom

% Number of ablated unit cells
nUCAblated = (ablationVolume .* (ucRatio / sum(ucRatio)) ./ ucVolume) .* densityRatio;

% Number of ablated atoms
nAtomAblated = 0;
for i = 1 : numel(uc)
    nAtomAblated = nAtomAblated + (nUCAblated(i) * sum(uc(i).AMOUNT));
end

% Pre-allocate memory
nParticleAngle = zeros(1, nAngle);

% Compute initial angular particle distribution (eq. 5)
for iAngle = 1 : (nAngle - 1)
    nParticleAngle(iAngle) = ((4/3)*pi * radiusDelta^3) ...
        .* abs( cosd(angle(iAngle + 1)) - cosd(angle(iAngle)) ) ...
        .* cosd(angle(iAngle)).^cosPowerFit;
end

% Normalize and multiply by the total number of ablated atoms
nParticleAngle = (nParticleAngle .* nAtomAblated) ./ sum(nParticleAngle);

%--------------------------------------------------------------------------
%% Write settings into config file
if createConfigBool
    createConfigFile( currentPath, folder, fileName, commentString, time, ...
        radius, angle, velo, bgPressure );
end

%--------------------------------------------------------------------------
%% Preallocation
% Number of property fields per velocity bin
nFieldsVelo = numel(enumeration('Field'));

% Memory of double [Bytes]
double = 8;

% Available RAM [GBytes]
memoryLimit = 4;

% Calculate memory size of plume and background particle matrices combined
% memorySize = nVelo * nFieldsVelo * nRadius * nTime * (nSpecies + 1) ...
%                 * nAngle * double * 10^-9;
memorySize = nVelo * nFieldsVelo * nRadius * (nAtomUCNumel + 1) * double * 10^-9;

% If the required memory exceeds the available memory, abort program and
% prompt error message
if (memorySize > memoryLimit)
    error(['Memory limit exceeded. Lower spatial, temporal, or angular' ...
            'resolution to reduce memory size of particle matrices.'])
end

% Plasma plume particle matrix
%   (velocity, field selector, radius, species)
plasmaMatrix = zeros(nVelo, nFieldsVelo, nRadius, nAtomUCNumel);

% Background gas particle matrix
%   (velocity, field selector, radius)
bgMatrix = zeros(nVelo, nFieldsVelo, nRadius);

%% Calculations
for iAngle = 1 : 1 % nAngle
    %% Main program per angle
    
    % If the number of atoms at current angle is smaller than the minimum
    %   number of particles per bin
    if nParticleAngle(iAngle) < nParticleMin
        % Skip to next angle
        continue
    end
    
    %% Fill initial background matrix
    for iRadius = 1 : (nRadius - 1)
        % Calculate bin volumes and background particle density per bin
        binVolume = (4/3)*pi ...
                    * ( (radius(iRadius + 1))^3 ...
                    - (radius(iRadius))^3 )     ...
                    * ( cosd(angle(iAngle))     ...
                    - cosd(angle(iAngle + 1)) ); % Bin volume [m^3]

        % Insert number of background particles per radial bin into
        %   first velocity bin (v = 0)
        bgMatrix(1, Field.nParticles, iRadius) = bgDensity * binVolume;
    end
    
    % Calculate the initial particle velocity distribution
    nParticleVeloInit = initialVelocityDistribution( plotVelDisInit, ...
        saveVelDisInit, velo, veloDistributionWidth, nUCAblated, uc, ...
        ucRatio, energyLaser, heatTarget, absorption, ...
        nParticleAngle(iAngle), densityRatio );
    
    % Only plot and save initial particle velocity distribution for the
    % center of the plume
    if iAngle == 1
        plotVelDisInit = false;
        saveVelDitInit = false;
    end
    
    % Fill first radial bin of plume particle matrix with the initial
    %   number of particles per angle per velocity bin
    for species = 1 : nAtomUCNumel
        plasmaMatrix(:, Field.nParticles, 1, species) = ...
            nParticleVeloInit(:, species);
    end
    
    for iTime = 1 : nTime
        %% Main program per timestep
        
        for iRadius = 1 : (nRadius - 1)
            %% Main program per radial bin
            
            % Loop though velocity bins from highest to lowest
            for iVelo = nVelo : -1 : 1
                %% Main program per velocity bin
                
                % If no radius bin is traveled this time step, the
                %   particles are assumed static
                if nRadiusBinTraveled(iVelo) == 0
                    % Skip velocity bin
                    continue
                end
                
                % Limit number of traveled radial bins to the number of
                %   radial bins between current bin and last bin
                if nRadiusBinTraveled(iVelo) > (nRadius - iRadius)
                    nRadiusVelo = nRadius - iRadius;
                else
                    nRadiusVelo = nRadiusBinTraveled(iVelo);
                end
                
                % Loop through species in plasma
                for species = 1 : nAtomUCNumel
                    %% Main program per species
                    
                    % If the number of particles in current bin is larger
                    %   than the minimum number of particles
                    if (plasmaMatrix(iVelo, Field.nParticles, iRadius, ...
                            species) > nParticleMin)
                        
                        % Loop through distance traveled in steps radiusDelta
                        for iRadiusDelta = 1 : nRadiusVelo
                            % Calculate collision probability per radius bin


                        end % For radiusDelta
                        
                    end
                    
                end % For species in plume
                
            end % For velocity
            
        end % For radius
        
    end % For time
    
end % For angle
