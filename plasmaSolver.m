%% PLASMASOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            Author: Sam Borkent
%
% Based on paper by: Tom Wijnands
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Naming conventions:
%
%   bg      : Background
%   bin     : Computational bin with equal position and velocity
%   n       : Number of particles/bins
%   velo    : Velocity
%   uc      : Unit cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% Do not change!

%--------------------------------------------------------------------------
% Clear workspace
%--------------------------------------------------------------------------
clc, clear, close all

%--------------------------------------------------------------------------
% Path settings
%--------------------------------------------------------------------------

% Get current path
currentPath = pwd;

% Add function and class folders
addpath([currentPath '/functions/']);
addpath([currentPath '/classes/']);

%--------------------------------------------------------------------------
% Create class instances holding constants and material properties
%--------------------------------------------------------------------------
C   = PhysicalConstants;    % Physical constants
PT  = PeriodicTable;        % Atom properties
UC  = UnitCells;            % Material properties

%% User settings
% All user input goes here !!!

%--------------------------------------------------------------------------
% Preferences (true / false)
%--------------------------------------------------------------------------

% Enable / disable debug mode
%   * Used to calculate conservation of particles and check for negative
%       number of particles.
debugBool = false;

% Create a configuration file containing simulation settings
createConfigBool = false;

% Plot the initial velocity distribution
plotVelDisInit = true;

% Save the initial velocity distribution plot
saveVelDisInit = false;

%--------------------------------------------------------------------------
% File settings
%--------------------------------------------------------------------------

% Save location
folder      = 'results';        % Output folder
fileName    = 'config';         % File name of configuration file

% Add comment to output file
commentString = [ ' ' ...
    'Type your comment here.' ...
    ];

%--------------------------------------------------------------------------
% Computational restrictions
%--------------------------------------------------------------------------

% Minimal number of particles per bin
nMin = 1;

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 19E-6;     % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 4E-5;     % Radial step size [m]

% Velocity limits
veloMax         = 2.5E4; % Maximal initial velocity
veloStepsize    = 1;     % Velocity step size

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Unit cell
%   * See UnitCells
%   * Options: TiO2, SrTiO3, Li2O, LiO2, LiTi2O4, Li2TiO3, Li4Ti5O12
%   * Preface with UC.
%   * Any mixture of crystals is possible, eg. [TiO2 Li2O]
uc = UC.TiO2;

% Mixture ratio of species in target
%   * Ignore if target contains only one component
%   * Should be of same length as unit cell array
ucRatio = [1 3];

% Target density [kg / m^3]
%   * Set to 0 for single crystal target
targetDensity = 0;

% Absorption coefficient
%   Default: 0.6
absorption = 0.6;

%--------------------------------------------------------------------------
% Background gas parameters
%--------------------------------------------------------------------------

% Background gas pressure during depostion [Pa = 1E2 mbar]
bgPressure      = 0.1E2;

% Temperature of background gas [K]
%   Default: 300
bgTemperature   = 300;      

%--------------------------------------------------------------------------
% Laser parameters
%--------------------------------------------------------------------------

% Laser spot width [m]
spotWidth       = 1E-3;

% Laser spot height [m]
spotHeight      = 2.1E-3;

% Ablation depth into the target for a single pulse [m]
%   * For most accurate results calculate per target.
%   Default: 100E-9
ablationDepth   = 100E-9;

% Laser fluence [J / cm^2]
laserFluence    = 2.0;

% Heat dissipating into the target each pulse [J]
%   * Small for ceramic targets.
%   * Can be approximated by measuring the temperature increase of the
%       target after a fixed number of laser pulses, then use the specific
%       heat equation to calculate total heat dissipation, divide by the
%       number of pulses.
%   Default: 0
heatTarget      = 0;

%--------------------------------------------------------------------------
% Model parameters
%--------------------------------------------------------------------------

% Fit parameter of the cosines power in the angular plasma particle
%   distibution function
cosPowerFit = 16;

% Initial velocity distribution width
initVeloDisWidth = 1500;

%% Calculations
% Everything below is calculated automatically based on user input above.

%--------------------------------------------------------------------------
% Material calculations
%--------------------------------------------------------------------------

% Number of elements in target
nElements = numel(uc.ELEMENTS);

% % Set mixture ratio to 1 for one component target
% if numel(uc) == 1
% 	ucRatio = 1;
% end
% 
% % Check if the number of given materials in target is the same as the
% %   number of elements in the mixture ratio array
% if numel(uc) ~= numel(ucRatio)
%     error([ 'Number of elements of mixture ratio array must be the ' ...
%             'same as the unit cell array.' ]);
% end
% 
% % Initialize atom arrays
% atomUCTemp      = [];  % Temporary array for calculations
% nAtomUCTemp     = [];  % Array containing amount of each atom in the target
% ucVolume        = zeros(1, numel(uc)); % Array holding unit cell volumes
% 
% % Fill atom arrays based on the given unit cells
% for i = 1 : numel(uc)
%     atomUCTemp  = [atomUCTemp   uc(i).ELEMENTS];
%     nAtomUCTemp = [nAtomUCTemp  (uc(i).AMOUNT .* ucRatio(i))];
%     ucVolume(i) = uc(i).VOLUME;
% end
% 
% % Normalize atom amount array
% nAtomUCTemp = nAtomUCTemp ./ sum(ucRatio);
% nAtomUCTemp = nAtomUCTemp ./ min(nAtomUCTemp);
% 
% % Initialize atom number array
% elementNumberTemp = zeros(1, numel(atomUCTemp));
% 
% % Fill atom number array
% for i = 1 : numel(atomUCTemp)
%     elementNumberTemp(i) = atomUCTemp(i).NUMBER;
% end
% 
% % Get sorted unique atom array
% [elementNumber, elementIndex, ~] = unique(elementNumberTemp, 'sorted');
% 
% % Check if there is oxygen (Z = 8)
% if ~isempty(elementNumber(elementNumber == 8))
%     % Place oxygen at the end of the atom array
%     elementIndex(end + 1) = elementIndex(elementNumber == 8);
%     
%     % Remove from old array
%     elementIndex(elementNumber == 8) = [];
% end
% 
% % Insert unique atoms in array
% atomUC = atomUCTemp(elementIndex');
% 
% % Initialize amount of atom array
% nAtomUC = zeros(1, numel(elementNumber));
% 
% % Look through unique atoms
% for i = 1 : numel(elementNumber)    
%     % Insert the calculated ratio of each atom in atom amount array
%     nAtomUC(i) = sum(nAtomUCTemp(elementNumberTemp == atomUC(i).NUMBER));
% end

%--------------------------------------------------------------------------
% Deposition parameter calculations
%--------------------------------------------------------------------------

% Background gas density (ideal gas law)
bgDensity = bgPressure / (C.BOLTZMANN * bgTemperature);

% Ablation volume [m^3]
ablationVolume  = spotWidth * spotHeight * ablationDepth;

%--------------------------------------------------------------------------
% Laser parameter calculations
%--------------------------------------------------------------------------

% Laser energy [J]
energyLaser = (laserFluence * 10^4) * (spotWidth * spotHeight);

%--------------------------------------------------------------------------
% Atom calculations
%--------------------------------------------------------------------------

% Number of ablated unit cells
nUCAblated = ablationVolume / uc.VOLUME;

% Number of ablated atoms
nAtomAblated = nUCAblated * sum(uc.AMOUNT);

%% Axis creation

% Angular axis
angle   = angleMin : angleDelta : angleMax - angleDelta;
nAngle  = numel(angle);

% Temporal axis
time    = timeMin : timeDelta : timeMax - timeDelta;
nTime   = numel(time);

% Radial axis
radius  = radiusMin : radiusDelta : radiusMax - radiusDelta;
nRadius = numel(radius);

% Velocity array
veloDelta   = (radiusDelta / timeDelta) / veloStepsize;
velo        = 0 : veloDelta : veloMax - veloDelta;
nVelo       = numel(velo);

%% Pre-allocations

% Plasma particle matrix
plasmaMatrix = zeros(nVelo, nRadius);
oxideMatrix = plasmaMatrix;

% plasmaNCol = plasmaMatrix;

% Plasma particle update matrices
plasmaSub = zeros(nVelo, nRadius);  % Removing particles
plasmaAdd = plasmaSub;              % Adding particles

oxideAdd = plasmaSub;              % Adding oxidized particles

% Background particle matrix
bgMatrix = zeros(nVelo, nRadius);

% Background particle update matrices
bgSub = bgMatrix;   % Removing particles
bgAdd = bgMatrix;   % Adding particles

%% Pre-calculations

% Number of radial bins traveled in one time step per velocity bin
nRadiusTraveled = round( (velo .* timeDelta) ./ radiusDelta );

% Index of first velocity high enough to traverse one radial bin per
%   timestep
iFirstVelo = find(nRadiusTraveled, 1);

% Angular plasma particle distribution
nParticleAngle = angularDistribution( angle,        ...
                                      radiusDelta,  ...
                                      cosPowerFit,  ...
                                      nAtomAblated );

%% Initialize figures

% 1D plasma propagation figure
figPropagation1D = figure('Name', '1D plasma particle propagation');
hold on;

% 1D background propagation figure
figBGPropagation1D = figure('Name', '1D background particle propagation');
hold on;

% 1D plasma propagation figure
figOxidePropagation1D = figure('Name', '1D oxidized plasma particle propagation');
hold on;

%% Main program

% Temporary values for testing
sigma   = pi * ( 2*PT.O.RADIUS + PT.Ti.RADIUS )^2;
mass    = PT.Ti.MASS;
massBG  = 2*PT.O.MASS;
energyOxidation = 9.7795 * C.EV;

% Kinetic energy of the plasma particle
energyPlasma = (0.5 * mass) .* (velo.^2);
energyBG = (0.5 * massBG) .* (velo.^2);

% If the oxidation energy is smaller than the kinetic energy of the plasma
%   particle, the particle still has energy to spare for oxidation
iFirstOxidation = find(energyPlasma > energyOxidation, 1);
                                  
for iAngle = 1 : 1 % nAngle
%% Calculations per angle

% Skip angle if there are not enough particles present
if nParticleAngle(iAngle) < nMin
    continue
end

% Pre-allocate background gas particle matrix and fill matrix based on
%   calculated density, also retrieve bin volumes per radial bin
[bgMatrix, binVolume] = fillBGMatrix( bgDensity,    ...
                                      nVelo,        ...
                                      radius,       ...
                                      angle,        ...
                                      iAngle,       ...
                                      radiusDelta );
                                  
% Total number of background particles at this angle
nBGTotal = sum( bgMatrix(1, :) );

% Total number of plasma particles at this angle
nPlasmaTotal = nParticleAngle(iAngle) / sum(uc.AMOUNT);
                                  
% Calculate the initial particle velocity distribution
[nParticleVeloInit, veloPlasma] = initialVelocityDistribution(          ...
                                               plotVelDisInit,          ...
                                               saveVelDisInit,          ...
                                               velo,                    ...
                                               initVeloDisWidth,        ...
                                               nUCAblated,              ...
                                               uc,                      ...
                                               1,                       ...
                                               energyLaser,             ...
                                               heatTarget,              ...
                                               absorption,              ...
                                               nParticleAngle(iAngle),  ...
                                               1 );

% Reset plasma matrix
if iAngle > 1
    plasmaMatrix = plasmaMatrix.*0;
end
                                           
% Fill into the plasma matrix
plasmaMatrix(:, 1) = nParticleVeloInit(:, 1);

for iTime = 1 : nTime
%% Calculations per time step

%--------------------------------------------------------------------------
% Reset update matrices for t+1
%--------------------------------------------------------------------------

plasmaSub = plasmaSub.*0;
plasmaAdd = plasmaSub;
bgSub = bgMatrix.*0;
bgAdd = bgSub;

oxideAdd = oxideAdd.*0;

for iRadius = (nRadius - 1) : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

for iVelo = nVelo : -1 : iFirstOxidation
%% Calculations per oxidizing plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Only loop through bins which have high velocity, so oxidation can occur

%--------------------------------------------------------------------------
% Get plasma particle variables
%--------------------------------------------------------------------------

% Number of plasma particles in current bin
nPlasma = plasmaMatrix(iVelo, iRadius);

% Skip bins with number of particles below threshold
if nPlasma < nMin
    continue
end

% Temporary variable holding remaining number of traveling particles after
%   subsequent collisions with the background
nPlasmaTemp = nPlasma;

%--------------------------------------------------------------------------
% Set number of expected traveled radial bins if no collisions would occur
%--------------------------------------------------------------------------

% Restrict number of traveled bins to the total number of
%   radial bins
if ( iRadius + nRadiusTraveled(iVelo) ) > nRadius
    nRadiusDelta = nRadius - iRadius;
else
    nRadiusDelta = nRadiusTraveled(iVelo);
end

for jRadius = 0 : (nRadiusDelta - 1)
%% Calculations per traversed radial bin
% Loop from current bin to one bin before final bin
% Collisions occur first in the bin closest to starting bin

% Break loop if no more particles are left to collide
if nPlasmaTemp < nMin
    break
end

% Index of the current radius: starting bin + traveled bin
thisRadius = iRadius + jRadius;

for iVeloBG = 1 : (iVelo - 1)
%% Calculations per background velocity bin
% Only include background velocities smaller than the plasma velocity
% From slowest to highest as collisions with slow moving background
%   particles are more probable

% Break loop if no more particles are left to collide
if nPlasmaTemp < nMin
    break
end

%--------------------------------------------------------------------------
% Get number of background particles
%--------------------------------------------------------------------------

% Number of background particals in this bin
%   Original number minus the already collided particles
nBG = bgMatrix(iVeloBG, thisRadius) + bgSub(iVeloBG, thisRadius);

% Skip background velocity if the number of background particles is below
%   the threshold
if nBG < nMin
    continue
end

%--------------------------------------------------------------------------
% Calculate velocity weight factors
%--------------------------------------------------------------------------

% Relative velocity weight factors
veloWeight = (velo(iVelo) - velo(iVeloBG)) ...
          ./ (velo(iVelo) + velo(iVeloBG));

%--------------------------------------------------------------------------
% Calculate number of collisions
%--------------------------------------------------------------------------

% Background particle density in current bin
bgDensity = nBG / binVolume(thisRadius);

% Collision rate
colRate = bgDensity * sigma * radiusDelta * veloWeight;

% Number of collided particles
nCol = colRate * nPlasmaTemp;

% Skip iteration if the number of collisions is smaller than the threshold
if nCol < nMin
    continue
end

% Set maximum number of collisions
if nPlasmaTemp > nBG
    % Number of background particles is limiting factor
    nColMax = nBG;
else
    % Number of plasma particles is limiting factor
    nColMax = nPlasmaTemp;
end

% Limit the number of collisions to the maximum number of collisions
if nCol > nColMax
    nCol = nColMax;
end

% Update number of non-collided plasma particles
nPlasmaTemp = nPlasmaTemp - nCol;

%--------------------------------------------------------------------------
% Calculate new velocities after collision
%--------------------------------------------------------------------------

% Total kinetic energy (sum of plasma and background particle energies)
energyKinetic = energyPlasma(iVelo) + energyBG(iVeloBG);

% Net energy (total kinetic energy minus oxidation energy)
energyNet = energyKinetic - energyOxidation;

% New velocity after oxidation
newVelo = sqrt( 2 * energyNet / (mass + massBG) );

% New plasma particle velocity index
iNewVelo = round( newVelo / veloDelta );
              
% Prevent index out-of-range error
if iNewVelo < 1
    iNewVelo = 1;
elseif iNewVelo > nVelo
    iNewVelo = nVelo;
end

%--------------------------------------------------------------------------
% Calculate new radial positions after collision
%--------------------------------------------------------------------------

% New radial index of collided plasma particles
nNewRadiusDelta = round( velo(iNewVelo) * (nRadiusDelta - jRadius) ...
                         * timeDelta / (nRadiusDelta * radiusDelta) );
                       
% Prevent index out-of-range error
if thisRadius + nNewRadiusDelta > nRadius
    nNewRadiusDelta = nRadius - thisRadius;
end

%--------------------------------------------------------------------------
% Update matrices
%--------------------------------------------------------------------------

% Remove collided particles from starting position
plasmaSub(iVelo, iRadius) = ...
    plasmaSub(iVelo, iRadius) - nCol;
bgSub(iVeloBG, thisRadius) = ...
    bgSub(iVeloBG, thisRadius) - nCol;

% Add collided particles to new position and new velocity after oxidation
oxideAdd(iNewVelo, thisRadius + nNewRadiusDelta) = ...
    oxideAdd(iNewVelo, thisRadius + nNewRadiusDelta) + nCol;

end % Background velocity loop

end % Traversed radial bin loop

end % Oxidizing plasma velocity loop

for iVelo = (iFirstOxidation - 1) : -1 : iFirstVelo
%% Calculations per non-oxidizing plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Only loop through velocities that traverse at least one radial bin per
%   time step
% Only loop through bins which have low velocity, so no oxidation can occur

%--------------------------------------------------------------------------
% Get plasma particle variables
%--------------------------------------------------------------------------

% Number of plasma particles in current bin
nPlasma = plasmaMatrix(iVelo, iRadius);

% Skip bins with number of particles below threshold
if nPlasma < nMin
    continue
end

% Temporary variable holding remaining number of traveling particles after
%   subsequent collisions with the background
nPlasmaTemp = nPlasma;

% kPlasma = plasmaNCol(iVelo, iRadius);

%--------------------------------------------------------------------------
% Set number of expected traveled radial bins if no collisions would occur
%--------------------------------------------------------------------------

% Restrict number of traveled bins to the total number of
%   radial bins
if ( iRadius + nRadiusTraveled(iVelo) ) > nRadius
    nRadiusDelta = nRadius - iRadius;
else
    nRadiusDelta = nRadiusTraveled(iVelo);
end

for jRadius = 0 : (nRadiusDelta - 1)
%% Calculations per traversed radial bin
% Loop from current bin to one bin before final bin
% Collisions occur first in the bin closest to starting bin

% Break loop if no more particles are left to collide
if nPlasmaTemp < nMin
    break
end

% Index of the current radius: starting bin + traveled bin
thisRadius = iRadius + jRadius;

% Total number of particles in this radial bin
% thisNPlasma = sum( plasmaMatrix(:, thisRadius) );   % Plasma
% thisNBG = sum( bgMatrix(:, thisRadius) );           % Background

% [DEBUG]
if debugBool == true
    % Reset total number of collisions in this radial bin
    nColSum = 0;
end

for iVeloBG = 1 : (iVelo - 1)
%% Calculations per background velocity bin
% Only include background velocities smaller than the plasma velocity
% From slowest to highest as collisions with slow moving background
%   particles are more probable

% Break loop if no more particles are left to collide
if nPlasmaTemp < nMin
    break
end

%--------------------------------------------------------------------------
% Get number of background particles
%--------------------------------------------------------------------------

% Number of background particals in this bin
%   Original number minus the already collided particles
nBG = bgMatrix(iVeloBG, thisRadius) + bgSub(iVeloBG, thisRadius);

% Skip background velocity if the number of background particles is below
%   the threshold
if nBG < nMin
    continue
end

%--------------------------------------------------------------------------
% Calculate velocity weight factors
%--------------------------------------------------------------------------

% Relative velocity weight factors
veloWeight = (velo(iVelo) - velo(iVeloBG)) ./ (velo(iVelo) + velo(iVeloBG));

%--------------------------------------------------------------------------
% Calculate number of collisions
%--------------------------------------------------------------------------

% Background particle density in current bin
bgDensity = nBG / binVolume(thisRadius);

% Normalization factor for the collision rate:
%   Scale by the fraction of non-colliding plasma particles compared to the
%   total number of plasma particles in this radial bin
% colNorm = nPlasmaTemp / (nPlasmaTemp + thisNPlasma);

% Collision rate
% colRate = bgDensity * sigma * radiusDelta * colNorm * veloWeight;
colRate = bgDensity * sigma * radiusDelta * veloWeight;

% [DEBUG]
if debugBool == true
    if colRate > 1
        disp('Number of collisions greater than 1.');
    end
end

% Normalization factor for the number of collisions:
%   Scale by the proportion of background particles of current verlocity
%   compared to all background particles in current radial bin
% nColNorm = nBG / thisNBG;

% Number of collided particles
% nCol = colRate * nPlasmaTemp * nColNorm;
nCol = colRate * nPlasmaTemp;

% Skip iteration if the number of collisions is smaller than the threshold
if nCol < nMin
    continue
end

% Set maximum number of collisions
if nPlasmaTemp > nBG
    % Number of background particles is limiting factor
    nColMax = nBG;
else
    % Number of plasma particles is limiting factor
    nColMax = nPlasmaTemp;
end

% Limit the number of collisions to the maximum number of collisions
if nCol > nColMax
    nCol = nColMax;
end

% Update number of non-collided plasma particles
nPlasmaTemp = nPlasmaTemp - nCol;

%--------------------------------------------------------------------------
% Calculate new velocities after collision
%--------------------------------------------------------------------------

% New plasma particle velocity
iNewVelo = round( ( iVelo*(mass - massBG) + 2*massBG*iVeloBG ) ...
                  / (mass + massBG) );

% New background particle velocity
iNewVeloBG = round( (iVeloBG*(massBG - mass) + 2 * mass * iVelo) ...
                     / (mass + massBG) );
              
% Prevent index out-of-range error
if iNewVelo < 1
    iNewVelo = 1;
elseif iNewVelo > nVelo
    iNewVelo = nVelo;
end
if iNewVeloBG < 1
    iNewVeloBG = 1;
elseif iNewVeloBG > nVelo
    iNewVeloBG = nVelo;
end

%--------------------------------------------------------------------------
% Calculate new radial positions after collision
%--------------------------------------------------------------------------

% New radial index of collided plasma particles
nNewRadiusDelta = round( velo(iNewVelo) * (nRadiusDelta - jRadius) ...
                         * timeDelta / (nRadiusDelta * radiusDelta) );
                       
% New radial index of collided background particles
nNewRadiusDeltaBG = round( velo(iNewVeloBG) * (nRadiusDelta - jRadius) ...
                           * timeDelta / (nRadiusDelta * radiusDelta) );

% % New radial index of collided plasma particles
% nNewRadiusDelta = round( (velo(iVelo) * jRadius ...
%                          + velo(iNewVelo) * (nRadiusDelta - jRadius) ) ...
%                          * timeDelta / (nRadiusDelta * radiusDelta) );
%                        
% % New radial index of collided background particles
% nNewRadiusDeltaBG = round( (velo(iVeloBG) * jRadius ...
%                            + velo(iNewVeloBG) * (nRadiusDelta - jRadius) ) ...
%                            * timeDelta / (nRadiusDelta * radiusDelta) );
                       
% Prevent index out-of-range error
if thisRadius + nNewRadiusDelta > nRadius
    nNewRadiusDelta = nRadius - thisRadius;
end
% if iRadius + nNewRadiusDelta > nRadius
%     nNewRadiusDelta = nRadius - iRadius;
% end
if thisRadius + nNewRadiusDeltaBG > nRadius
    nNewRadiusDeltaBG = nRadius - thisRadius;
end

%--------------------------------------------------------------------------
% Update matrices
%--------------------------------------------------------------------------

% Remove collided particles from starting position
plasmaSub(iVelo, iRadius) = ...
    plasmaSub(iVelo, iRadius) - nCol;
bgSub(iVeloBG, thisRadius) = ...
    bgSub(iVeloBG, thisRadius) - nCol;

% Add collided particles to new position and new velocity after collision
plasmaAdd(iNewVelo, thisRadius + nNewRadiusDelta) = ...
    plasmaAdd(iNewVelo, thisRadius + nNewRadiusDelta) + nCol;
bgAdd(iNewVeloBG, thisRadius + nNewRadiusDeltaBG) = ...
    bgAdd(iNewVeloBG, thisRadius + nNewRadiusDeltaBG) + nCol;

% % Add collided particles to new position and new velocity after collision
% plasmaAdd(iNewVelo, iRadius + nNewRadiusDelta) = ...
%     plasmaAdd(iNewVelo, iRadius + nNewRadiusDelta) + nCol;
% bgAdd(iNewVeloBG, thisRadius + nNewRadiusDeltaBG) = ...
%     bgAdd(iNewVeloBG, thisRadius + nNewRadiusDeltaBG) + nCol;

% plasmaNCol(iNewVelo, iRadius + nNewRadiusDelta) = kPlasma + 1;

%--------------------------------------------------------------------------
% Debug per collision event
%--------------------------------------------------------------------------

if debugBool == 1
    % [DEBUG] Total number of collisions in this radial bin
    nColSum = nColSum + nCol;
    
    % If the total number of collisions in this radial bin is greater than
    %   the total number of background particles in this bin
    if (nColSum - thisNBG) > nMin
    disp(['Too many background collisions: ' num2str(nColSum - thisNBG, '%.2E') '( ' ...
          num2str(iTime) ' ' num2str(iRadius) ' ' num2str(iVelo) ' ' ...
          num2str(thisRadius) ' ' num2str(iVeloBG) ' )' ]);
    % If the total number of collisions in this radial bin is greater the
    %   number of plasma particles before any collisions
    elseif nColSum > nPlasma
        disp('Too many plasma collisions.');
    end
end

end % Background velocity loop

end % Traversed radial bin loop

end % Non-oxidizing plasma velocity loop

end % Radius loop

%--------------------------------------------------------------------------
% Update non-colliding particles
%--------------------------------------------------------------------------

% Remove previously moved collided particles
plasmaMatrix = plasmaMatrix + plasmaSub;
bgMatrix = bgMatrix + bgSub;

% Update non-collided particles
plasmaMatrix = updateMatrix( plasmaMatrix, nRadius, nVelo, ...
                             iFirstVelo, nMin, nRadiusTraveled );
bgMatrix = updateMatrix( bgMatrix, nRadius, nVelo, ...
                         iFirstVelo, nMin, nRadiusTraveled );
                     
oxideMatrix = updateMatrix( oxideMatrix, nRadius, nVelo, ...
                            iFirstVelo, nMin, nRadiusTraveled );

% Add back the previously moved collided particles
plasmaMatrix = plasmaMatrix + plasmaAdd;
bgMatrix = bgMatrix + bgAdd;

oxideMatrix = oxideMatrix + oxideAdd;

%--------------------------------------------------------------------------
% Plot 1D propagation
%--------------------------------------------------------------------------

% Only for first angle (center of the plume)
if iAngle == 1   
    % Only for specific times
    if (iTime == 11)  || (iTime ==  16) || (iTime == 31) || (iTime ==  46) || ...
       (iTime == 51)  || (iTime ==  61) || (iTime == 71)
%     if (iTime == 11)  || (iTime ==  21) || (iTime == 31) || (iTime ==  41) || ...
%        (iTime == 51)  || (iTime ==  61) || (iTime == 91) || (iTime == 121) || ...
%        (iTime == 151) || (iTime == 181)

        nOxideTotal = sum(sum( oxideMatrix ));

        % Calculate total number of particles per radial bin
        nPlasmaRadius = nParticlesPerRadius( radius, plasmaMatrix );
        nBGRadius = nParticlesPerRadius( radius, bgMatrix );
        
        nOxideRadius = nParticlesPerRadius( radius, oxideMatrix );
        
        % Smooth data
        nPlasmaRadius = smoothdata(nPlasmaRadius, 2, 'gaussian', 100);
        nBGRadius = smoothdata(nBGRadius, 2, 'gaussian', 50);
        
        nOxideRadius = smoothdata(nOxideRadius, 2, 'gaussian', 50);
                
        % Normalization factor to conserve number of particles
        nPlasmaNorm = nPlasmaTotal / sum(nPlasmaRadius);
        nBGNorm = nBGTotal / sum(nBGRadius);
        
        nOxideNorm = nOxideTotal / sum(nOxideRadius);
        
        % Normalize
        nPlasmaRadius = nPlasmaRadius .* nPlasmaNorm;
        nBGRadius = nBGRadius .* nBGNorm;
        
        nOxideRadius = nOxideRadius .* nOxideNorm;
        
        % Fit a normalized Gaussian curve to the number of particles
%         nPlasmaFit = fitGaussian( radius, nPlasmaRadius, nMin );
        
        % Multiply by the total number of particles at this angle
%         nPlasmaRadius = nPlasmaFit .* nParticleAngle(iAngle);

        % Plot 1D plasma propagation
        figure(figPropagation1D);
        plot( radius, nPlasmaRadius, ...
              'LineWidth', 2, ...
              'DisplayName', [num2str(time(iTime), 3) ' s' ] );
          
        % Plot 1D background propagation
        figure(figBGPropagation1D);
        plot( radius, nBGRadius ./ binVolume, ...
              'LineWidth', 2, ...
              'DisplayName', [num2str(time(iTime), 3) ' s' ] );
          
        % Plot 1D oxide propagation
        figure(figOxidePropagation1D);
        plot( radius, nOxideRadius, ...
              'LineWidth', 2, ...
              'DisplayName', [num2str(time(iTime), 3) ' s' ] );
        
    end
end

%--------------------------------------------------------------------------
% Debug per time step
%--------------------------------------------------------------------------

if debugBool == true
    % Calculate total number of particles for background and plasma
    nBGTotalNew = sum(sum( bgMatrix ));
    nPlasmaTotalNew = sum(sum( plasmaMatrix ));

    % Check conservation of background particles
    if (nBGTotal - nBGTotalNew) > nMin
        disp(['Background particles not conserved ' ...
              num2str(nBGTotal - nBGTotalNew, '%.3E') ...
              ' particles lost.']);
    elseif (nBGTotal - nBGTotalNew) < -nMin
        disp(['Background particles not conserved ' ...
              num2str(nBGTotalNew - nBGTotal, '%.3E') ...
              ' particles gained.']);
    end

    % Check conservation of plasma particles
    if (nPlasmaTotal - nPlasmaTotalNew) > nMin
        disp(['Plasma particles not conserved ' ...
              num2str(nParticleAngle(iAngle) - nPlasmaTotalNew, '%.3E') ...
              ' particles lost.']);
    elseif (nPlasmaTotal - nPlasmaTotalNew) < -nMin
        disp(['Plasma particles not conserved ' ...
              num2str(nPlasmaTotalNew - nParticleAngle(iAngle), '%.3E') ...
              ' particles gained.']);
    end
end

end % Time loop

end % Angle loop

%% Plot settings

% 1D plasma propagation
figure(figPropagation1D);
legend;
xlim([0 0.05]);
% ylim([0 5E12]);
title('1D propagation of Ti from TiO_2 target in 0.02 mbar O_2 background');
xlabel('Target distance [m]');
ylabel('Number of particles per radial bin');

% 1D background propagation
figure(figBGPropagation1D);
legend;
xlim([0 0.05]);
% ylim([0 1E21]);
title('1D propagation of background O_2 gas particles');
xlabel('Target distance [m]');
ylabel('Particle density [m^-^3]');

% 1D plasma propagation
figure(figOxidePropagation1D);
legend;
xlim([0 0.05]);
% ylim([0 5E12]);
title('1D propagation of TiO2 oxidized in 0.02 mbar O_2 background');
xlabel('Target distance [m]');
ylabel('Number of particles per radial bin');
