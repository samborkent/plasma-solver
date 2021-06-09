%% PLASMASOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Author: Sam Borkent
%
%   Based on paper by: Tom Wijnands
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Naming conventions:
%
%   bg      : Background
%   bin     : Computational bin with similar position, velocity, etc.
%   n       : Number of particles/bins
%   velo    : Velocity
%   uc      : Unit cell
%
% Species indexing:
%
%   No oxygen
%   1       : Background gas molecule
%   2       : 1st target component
%   :       
%   N       : Nth target component
%
%   Oxygen gas, no oxygen in target:
%   1       : Background gas molecule
%   2       : 1st target component
%   3       : 2nd target component
%   :
%   N       : Nth target component
%   N+1     : 1st target component, 1st oxide
%   :
%   N+O     : 1st target component, Oth oxide (given by nOxidePerElement)
%   N+O+1   : 2nd target component, 1st oxide
%   N+2O    : 2nd target component, Oth oxide
%   :
%   N+NO    : Nth target component, Oth oxide
%
%   Oxygen gas and oxygen in target:
%   1       : Background gas molecule
%   2       : Oxygen present in target
%   3       : 1st non-oxygen target component
%   4       : 2nd non-oxygen target component
%   :
%   N+1     : Nth non-oxygen target component
%   :
%   etc.      Continues same as situation with oxgen gas and no oxygen in target
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% Do not change!

%-------------------------------------------------------------------------------
% Clear workspace
%-------------------------------------------------------------------------------
clc, clear, close all

%-------------------------------------------------------------------------------
% Path settings
%-------------------------------------------------------------------------------

% Get current path
currentPath = pwd;

% Add function and class folders
addpath([currentPath '/functions/']);
addpath([currentPath '/classes/']);

%-------------------------------------------------------------------------------
% Create class instances holding constants and material properties
%-------------------------------------------------------------------------------

C   = PhysicalConstants;    % Physical constants
PT  = PeriodicTable;        % Atom properties
UC  = UnitCells;            % Material properties

%% User settings
% All user input goes here !!!

%-------------------------------------------------------------------------------
% File settings
%-------------------------------------------------------------------------------

% Output folder
folder      = 'results';

% File name of configuration file
fileName    = 'config';

% Add comment to configuration file
commentString = [ ' ' ...
    'Type your comment here.' ...
    ];

%-------------------------------------------------------------------------------
% Preferences (true / false)
%-------------------------------------------------------------------------------

% Enable / disable debug mode
%   * Checks conservation of number of particles
%   * Checks for negative number of particles
%   * Reduces performance
debugBool = true;

% Create a configuration file containing simulation settings
createConfigBool = false;

% Plot the initial velocity distribution
plotInitVeloDistBool = true;

% Plot density plot for background gas particles instead of number of particles
plotDensityBool = false;

% Plot figures in log scale
plotLogBool = false;

% Plot all angles (true) or just the centre of the plasma (false)
%   * Calculation for all angles is supported, but 2D plots have not been
%       implemented yet.
plot2DBool = false;

%-------------------------------------------------------------------------------
% Computational restrictions
%-------------------------------------------------------------------------------

% Minimal number of particles per bin
%   * Best way to improve performance
%   * Increases quantization error and visible noise in results
nMin    = 1E6;     % Plasma threshold
nMinBG  = 1E3;     % Background gas threshold

% Maximum number of collisions
%   * Currently only non-collided (k=0) and collided (k=1) results get plotted
%   * Separating the plots into more different number of collisions is
%       supported, but not implemented at this moment.
kMax = 1;

% Allow negative velocities
%   * Huge performance impact
%   * Should be optimized to improve performance
negativeVeloBool = false;

% Switch velocity smoothing for non-collided particles on or off (true / false)
%   * Signifiacnt impact on performance
%   * Reduces quantization error
%   * Conserves number of particles
veloSmoothNonColBool = false;

% Switch velocity smoothing for collided particles on or off (true / false)
%   * Less visible effect on results
%   * Huge impact on performance
veloSmoothColBool = false;

% Width of the Gaussian smoothing function for smoothing particle velocities
%   * Default: 3 (corresponds to ~400 m/s with default resolution)
smoothWidth = 3;

% Bound particles to the last radial bin, so the total number of particles is
%   conserved, and particles don't dissappear after passing the last radial bin
%   * Has a negative inpact on performance
%   * Only enable to check the conservation of particles, other than that it has
%       no physical advantage
keepParticleBool = false;

%-------------------------------------------------------------------------------
% Dimensional limits
%-------------------------------------------------------------------------------

% Indices of times to plot
plotTimes = [6 16 31 51 81 101];

% Temporal limits [s]
timeMax     = 10E-6;     % Simulation run time ( Default: 5-20E-6 )
timeDelta   = 1E-7;     % Time step duration  ( Default: 1E-7    )

% Angular limits [deg]
angleMax    = 90;       % End angle         ( Default: 90 )
angleDelta  = 3;        % Angular step size ( Default: 3  )

% Radial limits [m]
radiusMax   = 0.05;     % Target-substrate distance ( Default: 0.05 )
radiusDelta = 6E-5;     % Radial step size          ( Default: 6E-5 )

%-------------------------------------------------------------------------------
% Material parameters
%-------------------------------------------------------------------------------

% Unit cell
%   * Preface with UC.
%   * Any mixture of crystals is possible, eg. [UC.TiO2 UC.Li2O]
%   * To add new materials or to view material propertie, see UnitCells.m
%   Options: TiO2, SrTiO3, Li2O, LiO2, LiTi2O4, Li2TiO3, Li4Ti5O12
% uc = [UC.TiO2 UC.Li2O];
uc = UC.TiO2;

% Mixture ratio of species in target
%   * Ignore if target contains only one component
%   * Should be of same length as unit cell array
ucRatio = [3 1];

% Target density [kg / m^3]
%   * If set to 0, a perfect single crystal target is assumed
%   * For most accurate results measure per target
%   * Default: 0
targetDensity = 0;

% Ratio of laser energy that is absorbed by the target
%   * Determine per target for most accurate results
%   * Default: 0.6
% absorbedRatio = 0.77;
absorbedRatio = 0.6;

% The number of oxides that can form per element
%   * Example: From a TiO2 target TiO and TiO2 can form during propagation,
%       so nOxidesPerElement is 2.
%   * Default: 2
nOxidePerElement = 0;

% Heat dissipating into the target per laser pulse [J]
%   * Small for ceramic targets.
%   * Can be approximated using specific heat equation if material
%       properties are known
%   * Default: 0
heatTarget      = 0;

% Excitation energy per atom in unit cell [J]
%   * Default: 0
energyExcitation = 0;

%-------------------------------------------------------------------------------
% Background gas parameters
%-------------------------------------------------------------------------------

% Background gas molecule
bg = [PT.O PT.O];

% Background gas molecule mass [kg]
massBG = sum([bg.MASS]);

% Background gas pressure during depostion [Pa = 1E2 mbar]
%   * Default: 0.01-0.2E2;
bgPressure      = 0.1E2;

% Temperature of background gas [K]
%   * Default: 300
bgTemperature   = 300;

%-------------------------------------------------------------------------------
% Laser parameters
%-------------------------------------------------------------------------------

% Laser spot width [m]
%   * Default: 1E-3
spotWidth       = 1E-3;

% Laser spot height [m]
%   * Default: 2E-3
% spotHeight      = 2.3E-3;
spotHeight      = 2.1E-3;

% Ablation depth into the target for a single pulse [m]
%   * For most accurate results calculate per target
%   * Default: 100E-9
ablationDepth   = 100E-9;

% Laser fluence [J / cm^2]
%   * Default: 2.0
laserFluence    = 2;

%-------------------------------------------------------------------------------
% Model parameters
%-------------------------------------------------------------------------------

% Fit parameter of the cosines power in the angular plasma particle distibution
%   function
%  * Default: 25
cosPowerFit = 16;

% Initial velocity distribution type
%   * Choices: 'normal', 'log-normal', 'Maxwell-Boltzmann'
%   * Default: 'Maxwell-Boltzmann'
distributionType = 'Maxwell-Boltzmann';

% Initial velocity distribution width
%   * Only relevant for normal or log-normal distribution
%   * Default: 1500
initVeloDistWidth = 1645;

%% Calculations
% Everything below is calculated automatically based on user input above.

%-------------------------------------------------------------------------------
% Timer
%-------------------------------------------------------------------------------

% Time the entire script
durationTotalTimer = tic;

% Initialize measured times
%   * To prevent createConfigFile throwing an error.
durationTotal = 0;

%-------------------------------------------------------------------------------
% Material calculations
%-------------------------------------------------------------------------------

% Throw error if number of background particle threshold is higher than
%   number of plasma particles threshold.
if nMinBG > (nMin + 1)
    error(['Threshold for number of background particles should be lower or ' ...
           'equal to the threshold for number of plasma particles.']);
end

% If target consists of only one unit cell
if numel(uc) == 1
    % Set mixture ratio to 1
	ucRatio = 1;
    
    % Single crystal density
    scDensity = uc.DENSITY;
    
    % Unit cell volume
    ucVolume = uc.VOLUME;
else
    % Normalize unit cell ratio
    ucRatio = ucRatio ./ sum(ucRatio);
    
    % Weighted average single crystal density
    scDensity = sum( [uc.DENSITY] .* ucRatio );
    
    % Weighted average unit cell volume
    ucVolume = sum( [uc.VOLUME] .* ucRatio );
end

% Check if the number of given materials in target is the same as the
%   number of elements in the mixture ratio array
if numel(uc) ~= numel(ucRatio)
    error([ 'Number of elements of mixture ratio array must be the ' ...
            'same as the unit cell array.' ]);
end

% If target is single crystal
if targetDensity == 0
    % Density ratio is 1
    densityRatio = 1;
else
    % Density ratio is given by ratio between measured and weighted
    %   single crystal density
    densiyRatio = targetDensity / scDenstity;
    
    % If density ratio somehow exceeds 1, throw error
    if densityRatio > 1
        error('Density ratio greater than 1.');
    end
end

% Get unique elements in target and their normalized amounts
[atomUC, nAtomUC] = getUniqueElements(uc, ucRatio);

% % Number of elements in target
nElements = numel(atomUC);

% Number of possible species
nSpecies = calculateNSpecies(atomUC, bg, nOxidePerElement);

% Get particle masses
mass = [sum([bg.MASS]) atomUC.MASS];

% Get metals present in target
metalUC = atomUC([atomUC.NUMBER] ~= 8);

% If a non-oxygen element is present in the target and there is oxygen
%   present in the background gas
if (numel(metalUC) >= 1) && sum([bg.NUMBER] == 8)
    % Add oxides to mass array
    %   * Assumes oxides of form MO, MO2, MO3, etc.
    %   * This is only correct for some material, if a material forms
    %       different types of oxides it has to be added manually (e.g. Li2O)
    for iOxide = 1 : nOxidePerElement
        mass = [ mass ([metalUC.MASS] + iOxide*PT.O.MASS) ];
    end
end

% Transpose mass array
mass = mass';

% Get radii of plasma species
colCS = [atomUC.RADIUS];

if (numel(metalUC) > 0) && sum([bg.NUMBER] == 8)
    % Add radii of oxides
    %   * Assumes oxides of form MO, MO2, MO3, etc.
    %   * This is only correct for some material, if a material forms
    %       different types of oxides it has to be added manually (e.g. Li2O)
    for iOxide = 1 : nOxidePerElement
        colCS = [ colCS ([metalUC.RADIUS] + iOxide*PT.O.RADIUS) ];
    end
end

% Calculate collisions cross-sections
colCS = (pi .* (sum([bg.RADIUS]) + colCS) .^2)';

%-------------------------------------------------------------------------------
% Background gas parameter calculations
%-------------------------------------------------------------------------------

% Background gas density (ideal gas law) [m^-3]
bgDensity = bgPressure / (C.BOLTZMANN * bgTemperature);

%-------------------------------------------------------------------------------
% Laser parameter calculations
%-------------------------------------------------------------------------------

% Ablation volume [m^3]
ablationVolume  = spotWidth * spotHeight * ablationDepth;

% Laser energy [J]
energyLaser = (laserFluence * 10^4) * (spotWidth * spotHeight);
% energyLaser = 0.0299;

%-------------------------------------------------------------------------------
% Atom calculations
%-------------------------------------------------------------------------------

% Number of ablated unit cells
nUCAblated = (ablationVolume / sum([uc.VOLUME] .* ucRatio)) ...
             .* ucRatio .* densityRatio;

% Number of ablated atoms
nAtomAblated = 0;
for i = 1 : numel(uc)
    nAtomAblated = nAtomAblated + nUCAblated(i) .* sum([uc(i).AMOUNT]);
end

%-------------------------------------------------------------------------------
% Axis creation
%-------------------------------------------------------------------------------

% Angular axis
angle   = 0 : angleDelta : angleMax - angleDelta;

% Set number of angles
if plot2DBool
    nAngle  = numel(angle);
else
    nAngle = 1;
end

% Temporal axis
time    = 0 : timeDelta : timeMax - timeDelta;
nTime   = numel(time);
executionTime = zeros(1, nTime);

% Radial axis
radius  = 0 : radiusDelta : radiusMax - radiusDelta;
nRadius = numel(radius);
iRadiusRange = 1 : nRadius;

% Velocity step size
veloDelta = round( radiusDelta / timeDelta );

% Initial velocity axis
veloInit = round(0 : veloDelta : 10E4);

%-------------------------------------------------------------------------------
% Angular distribution
%-------------------------------------------------------------------------------

% Angular plasma particle distribution
nPlasmaAngle = angularDistribution( angle,        ...
                                    angleDelta,   ...
                                    radiusDelta,  ...
                                    cosPowerFit,  ...
                                    nAtomAblated );
                                  
%-------------------------------------------------------------------------------
% Plot settings
%-------------------------------------------------------------------------------

% Indices for number of collisions
%   (1) : Non-collided
%   (2) : Collided
%   (3) : Sum

% Array holding line colors for each number of collisions
colorArray = ['b' 'r' 'k'];

% Array holding display name for each number of collisions
plotArray = {'Uncollided', 'Collided', 'Total'};

% Figure index
plotTimeIndex = 0;

%% DEBUG

if debugBool
    % Total number of loops
    nLoops = 0;
end

%% Main program

for iAngle = 1 : nAngle
%% Calculations per angle

% Skip angle if there are not enough particles present
if nPlasmaAngle(iAngle) < nMin
    continue
end

%-------------------------------------------------------------------------------
% Initial plasma particle velocity distribution
%-------------------------------------------------------------------------------

% Calculate the initial particle velocity distribution
nParticleVeloInit = initialVelocityDistribution( uc,                    ...
                                                 ucRatio,               ...
                                                 atomUC,                ...
                                                 nUCAblated,            ...
                                                 veloInit,              ...
                                                 initVeloDistWidth,     ...
                                                 absorbedRatio,         ...
                                                 energyLaser,           ...
                                                 heatTarget,            ...
                                                 energyExcitation,      ...
                                                 plotInitVeloDistBool,  ...
                                                 nMin,                  ...
                                                 nPlasmaAngle(iAngle),  ...
                                                 distributionType );
                                             
% Only plot initial velocity distribution for the first angle
if iAngle > 1
    plotInitVeloDistBool = false;
end

%-------------------------------------------------------------------------------
% Create velocity axis
%-------------------------------------------------------------------------------

% Maximum required velocity
%   * Skip all velocities where the number of particles is below the
%       threshold
iVeloMax = find( any(nParticleVeloInit > nMin, 1) , 1, 'last');

% Calculate the maximum achieved velocity after collision
iVeloMaxAlt = round( ( veloInit(1).*(mass(1) - max(mass(2:end))) ...
                       + 2*max(mass(2:end))*veloInit(iVeloMax) ) ...
                     ./ ( mass(1) + max(mass(2:end)) ) ./ veloDelta );

% If the maximum possible velocity after collision is higher than the previous
%   obtained value
if iVeloMaxAlt > iVeloMax
    % Increase the maximum velocity
    iVeloMax = iVeloMaxAlt;
end

% Velocity axis
if negativeVeloBool
    velo = -veloInit(iVeloMax) : veloDelta : veloInit(iVeloMax);
else
    velo = 0 : veloDelta : veloInit(iVeloMax);
end

% Number of velocity bins
nVelo = numel(velo);

% Zero velocity index
if negativeVeloBool
    iVeloZero = round(nVelo / 2);
else
    iVeloZero = 1;
end

% Velocity index range
iVeloRange  = 1 : nVelo;

% Check if there is a zero velocity
if ~any(velo == 0)
    error('Dimensional error. No zero velocity is present.');
end

%-------------------------------------------------------------------------------
% Pre-allocations
%-------------------------------------------------------------------------------

% Particle matrix 
particleMatrix = zeros(nSpecies, kMax+1, nVelo, nRadius);

% Matrix holding new positions of collided particles
%   * Resets every time step
collisionMatrix = zeros(nSpecies, kMax, nVelo, nRadius);

%-------------------------------------------------------------------------------
% Fill background
%-------------------------------------------------------------------------------

% Fill matrix based on calculated density, and retrieve bin volume per radial bin
[particleMatrix(1, 1, iVeloZero, :), binVolume] = fillBGMatrix( bgDensity,   ...
                                                                radius,      ...
                                                                radiusDelta, ...
                                                                angle,       ...
                                                                angleDelta,  ...
                                                                iAngle );

% Total number of background particles at this angle
nBGTotal = sum( particleMatrix(1, 1, iVeloZero, :) );

%-------------------------------------------------------------------------------
% Fill plasma
%-------------------------------------------------------------------------------

% Fill into the plasma matrix
particleMatrix(2:nElements+1, 1, iVeloZero:end, 1) = ...
    reshape(nParticleVeloInit(:, 1:iVeloMax), nElements, 1, [], 1);

% Total number of plasma particles at this angle
nPlasmaTotal = nPlasmaAngle(iAngle);

% Total number of particles at this angle
nParticleTotal = nBGTotal + nPlasmaTotal;

%-------------------------------------------------------------------------------
% Pre-calculations
%-------------------------------------------------------------------------------
                             
% Calculate collision rate matrix
%   * Density of background particles, times length of one radial bin,
%       time the collision cross-section
%   * Prevents having to calculate the collision rate millions of times,
%       instead calculate the collision rate for every position at once
%       each time step
collisionRate = ( 1 ./ reshape(binVolume, 1, 1, 1, nRadius) ) ...
                .* radiusDelta .* colCS;
            
% Matrix containing all possible velocity weights for collision calculation
%   * Prevents having to calculate many times, instead just look up value
%       from matrix
veloWeight = abs( (velo' - velo) ./ (abs(velo') + abs(velo)) );

% % Matrices holding indices of new velocity of particles after collision
% iNewVeloMatrix = zeros(nSpecies-1, nVelo, nVelo);   % Plasma particles
% iNewVeloBGMatrix = zeros(nSpecies-1, nVelo, nVelo); % Background gas particles
% 
% % Loop through plasma species
% for iSpecies = 2 : nSpecies
%     % Calculate new velocity of plasma species after collision with a
%     %   background gas particle
%     iNewVeloMatrix(iSpecies-1, :, :) = reshape( iVeloZero + round( ...
%         ( (mass(iSpecies) - mass(1)).*velo' + 2*mass(1).*velo ) ...
%         ./ (mass(iSpecies) + mass(1)) ./ veloDelta), 1, nVelo, nVelo );
%     
%     % Calculate new velocity of background gas particle after collision with
%     %   plasma species
%     iNewVeloBGMatrix(iSpecies-1, :, :) = reshape( iVeloZero + round( ...
%         ( 2*mass(iSpecies).*velo' + (mass(1) - mass(iSpecies)).*velo ) ...
%         ./ (mass(iSpecies) + mass(1)) ./ veloDelta), 1, nVelo, nVelo );
% end
% 
% % Prevent index out-of-bounds error
% iNewVeloMatrix(iNewVeloMatrix < 1) = 1;
% iNewVeloBGMatrix(iNewVeloBGMatrix < 1) = 1;
% iNewVeloMatrix(iNewVeloMatrix > nVelo) = nVelo;
% iNewVeloBGMatrix(iNewVeloBGMatrix > nVelo) = nVelo;

for iTime = 1 : nTime
%% Calculations per time step

% Reset collision matrix each time step
collisionMatrix = collisionMatrix.*0;

% If number of particle is not conserved
if ~keepParticleBool
    % Recalculate number of particles each time step
    nBGTotal = sum( particleMatrix(1, :, :, :), 'all' );
    nPlasmaTotal = sum( particleMatrix(2:end, :, :, :), 'all' );
end

%-------------------------------------------------------------------------------
% Velocity smoothing
%-------------------------------------------------------------------------------

% Skip smoothing for the first time step as all velocity bins are filled
if iTime > 1
    if veloSmoothNonColBool
        % Smooth velocities of non-collided particles
        particleMatrix(2:end, 1, :, :) = ...
            smoothdata( particleMatrix(2:end, 1, :, :), 3, 'gaussian', smoothWidth );
    end
    
    if veloSmoothColBool
        % Smooth velocities of collided particles
        particleMatrix(:, 2:end, [1:iVeloZero-1 iVeloZero+1:end], :) = ...
            smoothdata( particleMatrix(:, 2:end, [1:iVeloZero-1 iVeloZero+1:end], :), 3, 'gaussian', smoothWidth );
    end
end

%-------------------------------------------------------------------------------
% Remove particle below threshold and normalize matrices to conserve number
%   of particles
%-------------------------------------------------------------------------------

% Remove all background gas particles below threshold
%   * nMinBG has to be lower than nMin
particleMatrix(particleMatrix < nMinBG) = 0;

% Remove all plasma particles below threshold
tempMatrix = particleMatrix(2:end, :, :, :);
tempMatrix(tempMatrix < nMin) = 0;
particleMatrix(2:end, :, :, :) = tempMatrix;

% Normalize number of particles and add removed particles back to filled bins
particleMatrix(1, :, :, :) = ( particleMatrix(1, :, :, :) ...
                               ./ sum(particleMatrix(1, :, :, :), 'all') ) ...
                             .* nBGTotal;
particleMatrix(2:end, :, :, :) = ( particleMatrix(2:end, :, :, :) ...
                                   ./ sum(particleMatrix(2:end, :, :, :), 'all') ) ...
                                 .* nPlasmaTotal;

for iRadius = flip( iRadiusRange( any( particleMatrix(2:end, :, :, 1:end-1) > nMin, [1 2 3] ) ) )
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice (flip)
% Skip last bin as particles cannot move further (1:end-1)
% Only loop through radial bins where any of the bins inside have a number
%   of plasma particles above the threshold

for iVelo = flip( iVeloRange( any( particleMatrix(2:end, :, :, iRadius) > nMin, [1 2] ) ) )
%% Calculations per non-oxidizing plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles (flip)
% Only loop through velocity bins where any of the bins inside have a number
%   of plasma particles above the threshold

% Skip zero velocity as no distance is traversed
if iVelo == iVeloZero
    continue
end

% Save number of plasma particles before any collisions occured
nPlasmaTemp = particleMatrix(2:end, :, iVelo, iRadius);

%-------------------------------------------------------------------------------
% Determine number of traveled bins without collisions
%-------------------------------------------------------------------------------

% If velocity if positive
if iVelo > iVeloZero
    % If the number of forwards traveled radial bins exceeds the last bin
    if ( iRadius + (iVelo - iVeloZero) ) > nRadius
        % Limit to distance to last bin
        nRadiusDelta = nRadius - iRadius;
    else
        % Set number of traveled bins based on velocity
        nRadiusDelta = iVelo - iVeloZero;
    end
% If velocity if negative
elseif iVelo < iVeloZero
    % If the number of backwards traveled radial bins exceeds the first bin
    if  ( iRadius + (iVelo - iVeloZero) ) < 1
        % Limit to distance to first bin
        nRadiusDelta = iRadius - 1;
    else
        nRadiusDelta = abs(iVelo - iVeloZero);
    end
end

% If no distance is traveled skip to next bin
if nRadiusDelta == 0
    continue
end

for jRadius = 0 : (nRadiusDelta - 1)
%% Calculations per traversed radial bin
% Loop from current bin to one bin before final bin
% Collisions occur first in the bin closest to starting bin

% Break loop if no more particles are left to collide
if all( nPlasmaTemp < nMin, 'all' )
    break
end

% Index of the current radius: starting bin + traveled bin
if iVelo > iVeloZero
    thisRadius = iRadius + jRadius;
elseif iVelo < iVeloZero
    thisRadius = iRadius - jRadius;
end

% Find the filled background particle velocity bins at this radius
if iVelo > iVeloZero
    iVeloBGArray = iVeloRange(1:iVelo-1);
    iVeloBGArray = iVeloBGArray( any(particleMatrix(1, :, 1:iVelo-1, thisRadius) > nMinBG, 2) );
elseif iVelo < iVeloZero
    iVeloBGArray = iVeloRange(iVelo+1:end);
    iVeloBGArray = iVeloBGArray( any(particleMatrix(1, :, iVelo+1:end, thisRadius) > nMinBG, 2) );
end

for iVeloBG = iVeloBGArray
%% Calculations per background velocity bin
% Only include background velocities smaller than the plasma velocity
% From slowest to highest as collisions with slow moving background
%   particles are more probable

% Break loop if no more particles are left to collide
if all( nPlasmaTemp < nMin, 'all' )
    break
end

% Number of background particles
nBG = sum(particleMatrix(1, :, iVeloBG, thisRadius), 2);

% If number of background particles is below threshold
if nBG < nMinBG
    continue
end

% Number of collisions array
nCol = collisionRate(:, 1, 1, thisRadius) ...   % Pre-calculated collision rate
       .* nBG...                                % Number of background particles
       .* veloWeight(iVelo, iVeloBG) ...        % Difference in velocity weight
       .* nPlasmaTemp;                          % Number of plasma particles left

% Skip iteration if the number of collisions is smaller than the threshold
if all( nCol < nMinBG, 'all' )
    continue
end

% If the number of collisions is greater than the number of background particles
if sum( nCol, 'all' ) > nBG
    % Normalize number of collisions to number of background particles
    %   * Minus 1 to avoid faulty comparisons
    nCol = (nCol ./ sum(nCol, 'all')) .* (nBG - 1);
end

% Update number of non-collided plasma particles
nPlasmaTemp = nPlasmaTemp - nCol;

%-------------------------------------------------------------------------------
% Calculate new velocity after collision
%-------------------------------------------------------------------------------

% Calculate the new velocity of plasma particles after collision
iNewVelo = iVeloZero + round( ( velo(iVelo).*(mass(2:end) - mass(1)) ...
                                + 2*mass(1)*velo(iVeloBG) ) ...
                              ./ ( mass(2:end) + mass(1) ) ./ veloDelta );

% Calculate the new velocity background gas particles after collision
iNewVeloBG = iVeloZero + round( ( velo(iVeloBG).*(mass(1) - mass(2:end)) ...
                                  + 2*mass(2:end)*velo(iVelo) ) ...
                                ./ ( mass(1) + mass(2:end) ) ./ veloDelta );

% Prevent index out-of-range error
iNewVelo(iNewVelo < 1) = 1;
iNewVelo(iNewVelo > nVelo) = nVelo;
iNewVeloBG(iNewVeloBG < 1) = 1;
iNewVeloBG(iNewVeloBG > nVelo) = nVelo;

%-------------------------------------------------------------------------------
% Calculate new radial positions after collision
%-------------------------------------------------------------------------------

% New radial index of collided plasma particles
nNewRadiusDelta = round( velo(iNewVelo) ...
                         .* (nRadiusDelta - jRadius) ...
                         ./ (nRadiusDelta * veloDelta) );
                       
% New radial index of collided background particles
nNewRadiusDeltaBG = round( velo(iNewVeloBG) ...
                           .* (nRadiusDelta - jRadius) ...
                           ./ (nRadiusDelta * veloDelta) );

% Prevent index out-of-range error
nNewRadiusDelta(thisRadius + nNewRadiusDelta < 1) = ...
    1 - thisRadius;
nNewRadiusDeltaBG(thisRadius + nNewRadiusDeltaBG < 1) = ...
    1 - thisRadius;
nNewRadiusDelta(thisRadius + nNewRadiusDelta > nRadius) = ...
    nRadius - thisRadius;
nNewRadiusDeltaBG(thisRadius + nNewRadiusDeltaBG > nRadius) = ...
    nRadius - thisRadius;

%-------------------------------------------------------------------------------
% Update collided particles
%-------------------------------------------------------------------------------

% % Remove plasma particles from initial position before collision
particleMatrix(2:end, :, iVelo, iRadius) = ...
    particleMatrix(2:end, :, iVelo, iRadius) - nCol;

% Remove background particles from initial position before collision
%   * The total number of collided particles get removed from the
%       non-collided and collided background gas particles bins based on
%       the ratio between non-collided and collided particles.
particleMatrix(1, :, iVeloBG, thisRadius) = ...
    particleMatrix(1, :, iVeloBG, thisRadius) ...
    - ( particleMatrix(1, :, iVeloBG, thisRadius) ...
        ./ sum(particleMatrix(1, :, iVeloBG, thisRadius), 2) ) ...
      * sum(nCol, 'all');

% Sum number of collisions over number of collisions per particle
nColSum = sum(nCol, 2);

% Loop through plasma species
for iSpecies = 2 : nSpecies
    % Add plasma particles to new position after collision
    collisionMatrix(iSpecies, 1, iNewVelo(iSpecies-1), thisRadius + nNewRadiusDelta(iSpecies-1)) = ...
        collisionMatrix(iSpecies, 1, iNewVelo(iSpecies-1), thisRadius + nNewRadiusDelta(iSpecies-1)) ...
        + nColSum(iSpecies-1);

    % Add background particles to new position after collision
    collisionMatrix(1, 1, iNewVeloBG(iSpecies-1), thisRadius + nNewRadiusDeltaBG(iSpecies-1)) = ...
        collisionMatrix(1, 1, iNewVeloBG(iSpecies-1), thisRadius + nNewRadiusDeltaBG(iSpecies-1)) ...
        + nColSum(iSpecies-1);
end

if debugBool
    % Increment number of total loops
    nLoops = nLoops + 1;
end

end % Background velocity loop

end % Traveled distance loop

end % Velocity loop

end % Radius loop

%-------------------------------------------------------------------------------
% Update non-collided particles
%-------------------------------------------------------------------------------

% Propagate non-collided particles
particleMatrix = updateMatrix( particleMatrix, nVelo, iVeloZero, keepParticleBool );

% Add collided particles back to particle matrix
particleMatrix(:, 2, :, :) = particleMatrix(:, 2, :, :) + collisionMatrix;

% Only for specific times
if (iAngle == 1) && ( sum(iTime == plotTimes) == 1 )
    
    % Increment figure index
    plotTimeIndex = plotTimeIndex + 1;
    
    % Loop through plasma species
    for iSpecies = 1 : nElements+1
        % Initialize figure
        fig = figure(iSpecies);
        
        % Set figure name
        if iSpecies == 1
            fig.Name = 'Background propagation';
        else
            fig.Name = [atomUC(iSpecies-1).SYMBOL ' propagation'];
        end
        
        % Select subplot
        subplot(1, numel(plotTimes), plotTimeIndex);
        hold on;
        
        % Set title
        title([num2str(time(iTime), 3) ' s']);
        
        % Initialize sum of all collisions array
        nParticleRadiusSum = zeros(1, nRadius);
        
        % Loop through number of collisions
        for k = 1 : kMax+1
            % Calculate total number of particles per radial bin
            nParticleRadius = sum( squeeze( particleMatrix(iSpecies, k, :, :) ), 1 );
            
            % If non-collided and plasma species
            if k == 1 && iSpecies ~= 1
                % Total number of particles for this number of collisions
                nParticlePerK = sum( particleMatrix(iSpecies, k, :, :), 'all' );

                % Fit a normalized Gaussian curve to the number of particles
                nParticleRadius = fitGaussian( radius, nParticleRadius, nMin )';

                % Multiply by number of particles
                nParticleRadius = nParticleRadius .* nParticlePerK;
            else
                % Smooth data
                nParticleRadius = smoothdata( nParticleRadius, 2, 'gaussian', round(nRadius / 100) );
            end
            
            % Add to sum of all collisions array
            nParticleRadiusSum = nParticleRadiusSum + nParticleRadius;
            
            % If background gas species and density plot is enabled
            if (iSpecies == 1) && plotDensityBool
                plot( radius, nParticleRadius ./ binVolume, ...
                      'Color', colorArray(k), ...
                      'LineWidth', 2, ...
                      'DisplayName', plotArray{k} );
            else
                plot( radius, nParticleRadius, ...
                      'Color', colorArray(k), ...
                      'LineWidth', 2, ...
                      'DisplayName', plotArray{k} );
            end
        end % Number of collisions loop
       
        % If background gas species and density plot is enabled
        if (iSpecies == 1) && plotDensityBool
            plot( radius, nParticleRadiusSum ./ binVolume, ...
                  'Color', colorArray(end), ...
                  'LineWidth', 2, ...
                  'DisplayName', plotArray{end} );
        else
            plot( radius, nParticleRadiusSum, ...
                  'Color', colorArray(end), ...
                  'LineWidth', 2, ...
                  'DisplayName', plotArray{end} );
        end
        
        % Set log-scale if selected
        if plotLogBool
            set(gca, 'YScale', 'log');
        end
        
        % For the first subplot
        if plotTimeIndex == 1
            % Set y-limit as maximal value in plot
            ylimMax = max(nParticleRadiusSum);
        end
        
        % Set y-limit of plot
        if iSpecies == 1
            ylim([nMinBG ylimMax]);
        else
            ylim([nMin ylimMax]);
        end
        
        % Enable legend for last subplot
        if plotTimeIndex == numel(plotTimes)
            legend;
            legend('boxoff');
            legend('Location', 'northeast');
        end

        hold off;
        
    end % Species loop
end % Time if

% Show time
disp(['Simulation time: ' num2str(time(iTime)) ' s (Execution time: ' num2str(toc(durationTotalTimer)) ' s)']);

% Save execution times in array
executionTime(iTime) = toc(durationTotalTimer);

end % Time loop

end % Angle loop

% Get total excecution time [s]
durationTotal = toc(durationTotalTimer);

%% Create config file

if createConfigBool
    createConfigFile(currentPath, folder, fileName, commentString, ...          % File settings
        saveFiguresBool, plot2DBool, ...                                        % Preferences
        nMin, kMax, veloSmoothNonColBool, veloSmoothColBool, smoothWidth, keepParticleBool, ... % Computational restrictions
        time, angle, radius, velo, ...                                          % Dimensional limits
        uc, ucRatio, targetDensity, absorbedRatio, ...                          % Material parameters
        bg, bgPressure, bgTemperature, bgDensity, ...                           % Background gas parameters
        spotWidth, spotHeight, ablationDepth, ablationVolume, laserFluence, energyLaser, heatTarget, ... % Laser parameters
        cosPowerFit, initVeloDistWidth, ...                                      % Model parameters
        nPlasmaTotal, nBGTotal, ...                                             % Number of particles
        durationTotal );                                                        % Timers                          
end
