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
    'Li4Ti5O12 Target at 0.1 mbar' ...
    ];

% Add label to saved figures
figureLabel = 'LTO 0.1 mbar';

%-------------------------------------------------------------------------------
% Preferences
%-------------------------------------------------------------------------------

% Create a configuration file containing simulation settings
createConfigBool = true;

% Enable / disable debug mode
%   * Checks conservation of number of particles
%   * Checks for negative number of particles
debugBool = true;

%-------------------------------------------------------------------------------
% Plot settings
%-------------------------------------------------------------------------------

% Indices of times to plot
plotTimes = [5 15 30 50 80 100];

% Plot the initial velocity distribution
plotVeloDistBool = true;

% Smooth plot results
smoothPlotBool = true;

% Plot density plot for background gas particles instead of number of particles
plotDensityBool = true;

% Plot figures in log scale
plotLogBool = true;

% Plot all angles (true) or just the centre of the plasma (false)
%   * Calculation for all angles is supported, but 2D plots have not been
%       implemented yet.
%   * If enabled the script will be repeated for each angle.
plot2DBool = false;

% Save all open figures when simulation is finished
saveFiguresBool = false;

% Format to save figures in
%   * Ex.: '.fig', '.svg', '.png', etc.
saveFormat = '.fig';

%-------------------------------------------------------------------------------
% Computational restrictions
%-------------------------------------------------------------------------------
% To improve performance, you want to decrease the number of computed loops,
%   while retaining high quality results.
% There is a exponential relationship between the number of loops and the
%   excecution time: (numbers based on low end i3-5005U CPU)
%       Number of loops     Execution time
%       *           1E6     ~2 min
%       *           1E7     ~30 min
%       *           1E8     ~30 million years

% Enable / disable collisions
collisionBool = false;

% Minimal number of particles per bin
%   * Best way to improve performance
%   * Don't increase the values too much, otherwise all particles will be skipped
%   * Experiment to maximize performance / quality of results
%   * nMinBG has te be smaller than or equal to nMin
%   * Default: 1E6 (Limit seems to be ~1E8)
nMin    = 1E5;     % Plasma threshold
nMinBG  = 1E5;     % Background gas threshold

% Maximum number of collisions
%   * Currently only non-collided (k=0) and collided (k=1) results get plotted
%   * Separating the plots into more different number of collisions is
%       supported, but not implemented at this moment.
kMax = 1;

% Allow negative velocities
%   * Huge performance impact
%   * Should be optimized to improve performance
negativeVeloBool = true;

% Switch velocity smoothing for non-collided particles on or off (true / false)
%   * Signifiacnt impact on performance
%   * Reduces quantization error
%   * Conserves number of particles
%   * Best method for increasing the quality of results
%   * Turn off for quick testing
veloSmoothNonColBool = true;

% Switch velocity smoothing for collided particles on or off (true / false)
%   * Less visible effect on results
%   * Huge impact on performance
veloSmoothColBool = true;

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

% Temporal limits [s]
timeMax     = 10E-6;    % Simulation run time       ( Default: 5-20E-6 )
timeDelta   = 1E-7;     % Time step duration        ( Default: 100 ns  )

% Angular limits [deg]
angleMax    = 90;       % End angle                 ( Default: 90 )
angleDelta  = 3;        % Angular step size         ( Default: 3  )

% Radial limits [m]
%   * Decreasing radiusDelta increases the resolution of the model
%   * Has significant impact on performance
radiusMax   = 0.05;     % Target-substrate distance ( Default: 50 mm  )
radiusDelta = 10E-5;     % Radial step size          ( Default: 4-6E-5 )

%-------------------------------------------------------------------------------
% Material parameters
%-------------------------------------------------------------------------------

% Unit cell
%   * Preface with UC.
%   * Any mixture of crystals is possible, eg. [UC.TiO2 UC.Li2O]
%   * To add new materials or to view material propertie, see UnitCells.m
%   Options: TiO2, SrTiO3, Li2O, LiO2, LiTi2O4, Li2TiO3, Li4Ti5O12
uc = UC.SrTiO3;

% Mixture ratio of species in target
%   * Ignore if target contains only one component
%   * Should be of same length as unit cell array
ucRatio = [3 1];

% The number of oxides that can form per element
%   * Example: From a TiO2 target TiO and TiO2 can form during propagation,
%       so nOxidesPerElement is 2.
%   * Default: 2
nOxidePerElement = 0;

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

% Background gas pressure during depostion [Pa = 1E2 mbar]
%   * Default: 0.01-0.2E2;
bgPressure      = 0.02E2;

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
cosPowerFit = 25;

% Initial velocity distribution type
%   * Choices: 'normal', 'log-normal', 'Maxwell-Boltzmann'
%   * Default: 'Maxwell-Boltzmann'
distributionType = 'Maxwell-Boltzmann';

% Initial velocity distribution width
%   * Only relevant for normal or log-normal distribution
%   * Default: 1500
initVeloDistWidth = 1500;

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

% Total number of computed loops
nLoops = 0;

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

atomUC = flip( atomUC );
nAtomUC = flip( nAtomUC );

% Number of elements in target
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

%-------------------------------------------------------------------------------
% Atom calculations
%-------------------------------------------------------------------------------

% Number of ablated unit cells
nUCAblated = (ablationVolume / sum([uc.VOLUME] .* ucRatio)) ...
             .* ucRatio .* densityRatio;

% Number of ablated atoms
nAtomAblated = 0;
for jSpecies = 1 : numel(uc)
    nAtomAblated = nAtomAblated + nUCAblated(jSpecies) .* sum([uc(jSpecies).AMOUNT]);
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
time    = timeDelta : timeDelta : timeMax;
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
                                                 plotVeloDistBool,  ...
                                                 nMin,                  ...
                                                 nPlasmaAngle(iAngle),  ...
                                                 distributionType );
                                             
% Only plot initial velocity distribution for the first angle
if iAngle > 1
    plotVeloDistBool = false;
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

% Matrix holding all collided particles that need to be subtracted from
%   particle matrix
subtractMatrix = zeros(nSpecies, kMax+1, nVelo, nRadius);

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
for iVelo = 1 : iVeloMax-1
    particleMatrix(2:nElements+1, 1, iVeloZero+iVelo, iVelo) = nParticleVeloInit(:, iVelo+1);
end

% Total number of plasma particles at this angle
nPlasmaTotal = nPlasmaAngle(iAngle);

% Total number of particles at this angle
nParticleTotal = nBGTotal + nPlasmaTotal;

%-------------------------------------------------------------------------------
% Calculate initial plasma and background gas temperatures
%-------------------------------------------------------------------------------

[plasmaTempStart, bgTempStart] = averageTemperature( particleMatrix, velo, mass );

for iTime = 1 : nTime
%% Calculations per time step

% Reset collision matrix each time step
collisionMatrix = collisionMatrix.*0;
subtractMatrix = subtractMatrix.*0;

% If number of particle is not conserved
if ~keepParticleBool
    % Recalculate number of particles each time step
    nBGTotalNew = sum( particleMatrix(1, :, :, :), 'all' );
    nPlasmaTotalNew = sum( particleMatrix(2:end, :, :, :), 'all' );
end

%-------------------------------------------------------------------------------
% Velocity smoothing
%-------------------------------------------------------------------------------
% Avoids zero velocity and smooth positive and negative velocities seperately,
%   so no particles have their propagation direction reversed by smoothing,
%   and no particles become static due to smoothing

if veloSmoothNonColBool
    % Smooth velocities of non-collided particles
    particleMatrix(2:end, 1, 1:iVeloZero-1, :) = ...
        smoothdata( particleMatrix(2:end, 1, 1:iVeloZero-1, :), 3, 'gaussian', smoothWidth );
    particleMatrix(2:end, 1, iVeloZero+1:nVelo, :) = ...
        smoothdata( particleMatrix(2:end, 1, iVeloZero+1:nVelo, :), 3, 'gaussian', smoothWidth );
end

if veloSmoothColBool
    % Smooth velocities of collided particles
    particleMatrix(:, 2:end, 1:iVeloZero-1, :) = ...
        smoothdata( particleMatrix(:, 2:end, 1:iVeloZero-1, :), 3, 'gaussian', smoothWidth );
    particleMatrix(:, 2:end, iVeloZero+1:end, :) = ...
        smoothdata( particleMatrix(:, 2:end, iVeloZero+1:end, :), 3, 'gaussian', smoothWidth );
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
                             .* nBGTotalNew;
particleMatrix(2:end, :, :, :) = ( particleMatrix(2:end, :, :, :) ...
                                   ./ sum(particleMatrix(2:end, :, :, :), 'all') ) ...
                                 .* nPlasmaTotalNew;
                             
%-------------------------------------------------------------------------------
% Collision calculations
%-------------------------------------------------------------------------------
                             
% If collisions are enabled
if collisionBool
    % All collision calculation are performed here
    [particleMatrix, collisionMatrix, nLoops] = collisionCalculation( particleMatrix, ...
        nLoops, nMin, nMinBG, iRadiusRange, radiusDelta, nRadius, velo, veloDelta, nVelo, ...
        iVeloRange, iVeloZero, binVolume, colCS, mass, nSpecies );
end

%-------------------------------------------------------------------------------
% New stuff
%-------------------------------------------------------------------------------

% Number of loops counter
NLoops = 0;

% Number of possible types of collisions
N_C = 0;
for jSpecies = 1 : nSpecies-1
    N_C = N_C + (nSpecies - jSpecies);
end

% Initialize matrices
% v_AB = zeros(N_C, nRadius); % Weighted average velocity of particles A and B
% C_AB = zeros(N_C, nRadius); % Number of collisions between particles A and B
% C_AB = zeros(1, nVelo); % Number of collisions between particles A and B

% Radii of particles
particleRadius = [sum([bg.RADIUS]) atomUC.RADIUS];

% Loop through radii
for iRadius = flip( iRadiusRange( any( particleMatrix(2:end, :, :, :) > nMin, [1 2 3] ) ) )

% Direction
for dir = [1 -1]

% Loop through distance traveled
for iVelo = flip( (1 : iVeloMax-1).*dir )
% Type of collision index
index = 0;

% Loop through species
for iSpecies = nSpecies : -1 : 2
    
    % Skip loop if number of A particles is below threshold
    if all( particleMatrix(iSpecies, :, iVeloZero+iVelo, iRadius) < nMin, 2 )
        continue
    end

    % If the particle is moving forwards (towards substrate)
    if dir == 1
        % Set velocity range
        veloRange = 1:iVeloZero+iVelo-1;

        % Set radius range
        if iRadius+iVelo-1 > nRadius
            radiusRange = iRadius:nRadius;
        else
            radiusRange = iRadius:iRadius+iVelo-1;
        end
    % If the particle if moving backwards (towards target)
    elseif dir == -1
        % Set velocity range
        veloRange = iVeloZero+iVelo+1:nVelo;

        % Set radius range
        if iRadius+iVelo+1 < 1
            radiusRange = 1:iRadius;
        else
            radiusRange = iRadius+iVelo+1:iRadius;
        end
    end
    
    % Skip loop if all B particles are below threshold
    if all( particleMatrix( 1:iSpecies-1, :, veloRange, radiusRange ) < nMin, 'all' )
       continue
    end
    
    % Get all filled radius indices
    radiusRangeFull = radiusRange( sum(particleMatrix(1:iSpecies-1, :, veloRange, radiusRange), [1 2 3]) > nMin );

    % Get all filled velocity indices
    veloRangeFull = veloRange( sum(particleMatrix(1:iSpecies-1, :, veloRange, radiusRangeFull), [1 2 4]) > nMin );
    
    % Skip loop if ranges are empty
    if isempty(radiusRangeFull) || isempty(veloRangeFull)
        continue
    end
    
    % Relative velocities
    v_AB = abs( velo(iVeloZero+iVelo) - velo(veloRangeFull) );
    
    % Skip loop if relative velocity is smaller than the minimum non-zero velocity
    if all(round( v_AB ./ veloDelta ) < 1, 'all')
        continue
    end
    
    % Velocity of particle A after collision
    v_A = ( ( mass(iSpecies) - mass(1:iSpecies-1) ) * velo(iVeloZero+iVelo) ...
            + 2 * mass(1:iSpecies-1) .* velo(veloRangeFull) ) ...
          ./ ( mass(iSpecies) + mass(1:iSpecies-1) );

    % Velocity of particles B after collision
    v_B = ( ( mass(1:iSpecies-1) - mass(iSpecies) ) .* velo(veloRangeFull) ...
            + 2 * mass(iSpecies) .* velo(iVeloZero+iVelo) ) ...
          ./ ( mass(iSpecies) + mass(1:iSpecies-1) );
      
    % Velocity indices after collision
    v_A = iVeloZero + round( v_A ./ veloDelta );
    v_B = iVeloZero + round( v_B ./ veloDelta );  

    % Number of particles of type A and B
    N_A = sum( particleMatrix(iSpecies, :, iVeloZero+iVelo, iRadius), 2 );
    N_B = squeeze( sum( particleMatrix(1:iSpecies-1, :, veloRangeFull, radiusRangeFull), 2 ) );
    
    % Reshape so dimensions are correct
    N_B = reshape( N_B, numel(1:iSpecies-1), numel(veloRangeFull), numel(radiusRangeFull) );

    % Collisions cross-section of particles A and B
    sigma_AB = pi .* (particleRadius(iSpecies) + particleRadius(1:iSpecies-1))'.^2;
    
    % Number of collisions between particles A and B within one time step
    C_AB = ( N_A .* N_B ./ reshape(binVolume(radiusRangeFull), 1, 1, numel(radiusRangeFull) ) ) ...
        .* sigma_AB .* timeDelta .* v_AB;

    % Remove value if number of collisions is below threshold
    C_AB(C_AB < nMin) = 0;

    % Skip loop if all number of collisions have been set to zero
    if all(C_AB == 0, 'all')
        continue
    end

    % Limit the number of collisions to the number of A particles
    if sum(C_AB, 'all') > N_A
        C_AB = ( C_AB ./ sum(C_AB, 'all') ) .* N_A;
    end

    % Limit the number of collisions to the number of B particles
    C_AB( C_AB > N_B ) = N_B( C_AB > N_B );

    % If there are non-collided A particles
    if particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) > nMin
        
        % If there are enough non-collided A particles
        if particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) > sum(C_AB, 'all')
            
            % Remove non-collided A particles
            subtractMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) = ...
                subtractMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) + sum(C_AB, 'all');

        % If there are not enough non-collided A particles
        else
            % Remove remainder from collided A particles
            subtractMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) = ...
                subtractMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) ...
                + sum(C_AB, 'all') ...
                - particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius);
            
            % Remove all non-collided A particles
            subtractMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) = ...
                subtractMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) ...
                + particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius);
        end
    % If there are no non-collided A particles left
    else
        % Remove all collisions from collided A particles
         subtractMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) = ...
             subtractMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) + sum(C_AB, 'all');
    end
    
    % Temporary B particle matrix
    bTemp = squeeze( particleMatrix(1:iSpecies-1, 1, veloRangeFull, radiusRangeFull) );
    
    bTemp = reshape( bTemp, size(C_AB) );
    cols1 = C_AB .* 0;
    cols2 = C_AB .* 0;
    
    % If there are non-collided B particles left
    if any( particleMatrix(1:iSpecies-1, 1, veloRangeFull, radiusRangeFull) > nMin, 'all' )
        
        % If there are enough non-collided B particles
        if any( bTemp > C_AB, 'all' )
            
            cols1(bTemp > C_AB) = C_AB(bTemp > C_AB);
            cols1(bTemp < C_AB) = bTemp(bTemp < C_AB);
            cols2(bTemp < C_AB) = C_AB(bTemp < C_AB) - bTemp(bTemp < C_AB);
                     
            % Remove non-collided B particles
            subtractMatrix(1:iSpecies-1, 1, veloRangeFull, radiusRangeFull) = ...
                subtractMatrix(1:iSpecies-1, 1, veloRangeFull, radiusRangeFull) ...
                + reshape( cols1, numel(1:iSpecies-1), 1, numel(veloRangeFull), numel(radiusRangeFull) );
            
            % Remove non-collided B particles
            subtractMatrix(1:iSpecies-1, 2, veloRangeFull, radiusRangeFull) = ...
                subtractMatrix(1:iSpecies-1, 2, veloRangeFull, radiusRangeFull) ...
                + reshape( cols2, numel(1:iSpecies-1), 1, numel(veloRangeFull), numel(radiusRangeFull) );
                     
        % If there are not enough non-collided A particles
        else

            % Remove remainder from collided B particles
            subtractMatrix(1:iSpecies-1, 2, veloRangeFull, radiusRangeFull) = ...
                subtractMatrix(1:iSpecies-1, 2, veloRangeFull, radiusRangeFull) ...
                + reshape( C_AB - bTemp, numel(1:iSpecies-1), 1, numel(veloRangeFull), numel(radiusRangeFull) );
            
            % Remove all non-collided B particles
            subtractMatrix(1:iSpecies-1, 1, veloRangeFull, radiusRangeFull) = ...
                subtractMatrix(1:iSpecies-1, 1, veloRangeFull, radiusRangeFull) ...
                + particleMatrix(1:iSpecies-1, 1, veloRangeFull, radiusRangeFull);        
            
        end
    % If there are no non-collided B particles left
    else
        % Remove all collisions from collided B particles
        subtractMatrix(1:iSpecies-1, 2, veloRangeFull, radiusRangeFull) = ...
            subtractMatrix(1:iSpecies-1, 2, veloRangeFull, radiusRangeFull) ...
            + reshape( C_AB, numel(1:iSpecies-1), 1, numel(veloRangeFull), numel(radiusRangeFull) );
    end
    
%     % Scale the actually traveled distance by the difference in velocity
%     %   before and after collision
%     nRadiusTraveledA = round( 0.5.*(iVelo + ( velo(v_A) ./ velo(iVeloZero + iVelo) ) .* iVelo) );
    
    C_AB_Sum = sum(C_AB, 3);

    for jSpecies = 1 : iSpecies-1
        % Add collided A particles to collision matrix
        collisionMatrix(iSpecies, 1, v_A(jSpecies, :), iRadius+iVelo) = ...
            collisionMatrix(iSpecies, 1, v_A(1, :), iRadius+iVelo) ...
            + reshape( C_AB_Sum(jSpecies, :), 1, 1, numel(veloRangeFull), 1);
        
        % Add collided B particles to collision matrix
        collisionMatrix(jSpecies, 1, v_B(jSpecies, :), radiusRangeFull) = ...
            collisionMatrix(jSpecies, 1, v_B(jSpecies, :), radiusRangeFull) ...
            + reshape( C_AB(jSpecies, :, :), 1, 1, numel(veloRangeFull), numel(radiusRangeFull) );
        
        % Increment number of loops
        NLoops = NLoops + 1;
    end
    
% % Loop through species that are not current species
% for jSpecies = iSpecies-1 : -1 : 1                 
%     % Increment types of collisions index
%     index = index + 1;
% 
%     % If the particle is moving forwards (towards substrate)
%     if dir == 1
%         % Set velocity range
%         veloRange = 1:iVeloZero+iVelo-1;
% 
%         % Set radius range
%         if iRadius+iVelo-1 > nRadius
%             radiusRange = iRadius:nRadius;
%         else
%             radiusRange = iRadius:iRadius+iVelo-1;
%         end
%     % If the particle if moving backwards (towards target)
%     elseif dir == -1
%         % Set velocity range
%         veloRange = iVeloZero+iVelo+1:nVelo;
% 
%         % Set radius range
%         if iRadius+iVelo+1 < 1
%             radiusRange = 1:iRadius;
%         else
%             radiusRange = iRadius+iVelo+1:iRadius;
%         end
%     end
% 
%     % Skip loop if all particle are below threshold
%     if    all( particleMatrix(iSpecies, :, iVeloZero+iVelo, iRadius) < nMin, 2 ) ...
%        || all( particleMatrix(jSpecies, :, veloRange, radiusRange) < nMin, 'all' )
%        continue
%     end
% 
%     % Get all filled radius indices
%     radiusRangeFull = radiusRange( sum(particleMatrix(jSpecies, :, veloRange, radiusRange), [2 3]) > nMin );
% 
%     % Get all filled velocity indices
%     veloRangeFull = veloRange( sum(particleMatrix(jSpecies, :, veloRange, radiusRangeFull), [2 4]) > nMin );
%     
%     if isempty(radiusRangeFull) || isempty(veloRangeFull)
%         continue
%     end
%     
%     % Velocity of particle A after collision
%     v_A = ( ( mass(iSpecies) - mass(jSpecies) ) * velo(iVeloZero+iVelo) ...
%             + 2 * mass(jSpecies) .* velo(veloRangeFull) ) ...
%           ./ ( mass(iSpecies) + mass(jSpecies) );
% 
%     % Velocity of particle B after collision
%     v_B = ( ( mass(jSpecies) - mass(iSpecies) ) .* velo(veloRangeFull) ...
%             + 2 * mass(iSpecies) .* velo(iVeloZero+iVelo) ) ...
%           ./ ( mass(iSpecies) + mass(jSpecies) );
% 
%     % Velocity indices after collision
%     v_A = iVeloZero + round( v_A ./ veloDelta )';
%     v_B = iVeloZero + round( v_B ./ veloDelta )';
% 
%     % Number of particles of type A and B
%     N_A = sum( particleMatrix(iSpecies, :, iVeloZero+iVelo, iRadius), 2 );
%     N_B = squeeze( sum( particleMatrix(jSpecies, :, veloRangeFull, radiusRangeFull), 2 ) );
%     
%     % If the dimensions don't agree
%     if [numel(veloRangeFull) numel(radiusRangeFull)] ~= size(N_B)
%         % Transpose array
%         N_B = N_B';
%     end
% 
%     % Average velocity of particles A and B combined
%     v_AB = abs( velo(iVeloZero+iVelo) - velo(veloRangeFull) )';
%     
%     % Skip loop if it is smaller than the minimum non-zero velocity
%     if round( v_AB / veloDelta ) < 1
%         continue
%     end
% 
%     % Radius of particle B
%     r_A = atomUC(iSpecies-1).RADIUS;
% 
%     % Radius of particle A
%     if jSpecies == 1
%         r_B = bg.RADIUS;
%     else
%         r_B = atomUC(jSpecies-1).RADIUS;
%     end
% 
%     % Collisions cross-section of particles A and B
%     sigma_AB = pi * (r_A + r_B)^2;
% 
%     % Number of collisions between particles A and B within one time step
%     C_AB = ( (N_A / iVelo) .* N_B ./ binVolume(radiusRangeFull) ) .* sigma_AB .* timeDelta .* v_AB;
% 
%     % Total number of collisions
%     C_AB_Sum = sum(C_AB, 'all');
% 
%     % Remove value if number of collisions is below threshold
%     if any(C_AB < nMin, 'all')
%         C_AB(C_AB < nMin) = 0;
%     end
% 
%     % Skip loop if all number of collisions have been set to zero
%     if all(C_AB == 0, 'all')
%         continue
%     end
% 
%     % Normalize number of collisions
%     C_AB = (C_AB ./ sum(C_AB, 'all')) .* C_AB_Sum;
% 
%     % Limit the number of collisions to the number of A particles
%     if sum(C_AB, 'all') > N_A
%         C_AB = ( C_AB ./ sum(C_AB, 'all') ) .* N_A;
%     end
%     
%     % Limit the number of collisions to the number of B particles
%     if any( C_AB > N_B, 'all' )
%         C_AB( C_AB > N_B ) = N_B( C_AB > N_B );
%     end
%     
%     % If there are non-collided A particles
%     if particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) > nMin
%         
%         % If there are enough non-collided A particles
%         if particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) > sum(C_AB, 'all')
%             
%             % Remove non-collided A particles
%             particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) = ...
%                 particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) - sum(C_AB, 'all');
%             
%             if any(particleMatrix < -nMin, 'all')
%                 disp('There are non-collided A particles, and there are more then collisions.');
%                 disp([iTime iRadius iVeloZero+iVelo iSpecies jSpecies]);
%                 return
%             end
%         
%         % If there are not enough non-collided A particles
%         else
%             
%             % Remove remainder from collided A particles
%             particleMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) = ...
%                 particleMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) ...
%                 - ( sum(C_AB, 'all') - particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) );
%             
%             % Remove all non-collided A particles
%             particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) = 0;
%             
%             if any(particleMatrix < -nMin, 'all')
%                 disp('There are non-collided A particles, but there are more collisions.');
%                 disp([iTime iRadius iVeloZero+iVelo iSpecies jSpecies]);
%                 return
%             end
%             
%         end
%     % If there are no non-collided A particles left
%     else        
%         % Remove all collisions from collided A particles
%          particleMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) = ...
%             particleMatrix(iSpecies, 2, iVeloZero+iVelo, iRadius) - sum(C_AB, 'all');
%         
%         if any(particleMatrix < -nMin, 'all')
%             disp('There are no non-collided A particles.');
%             disp([iTime iRadius iVeloZero+iVelo iSpecies jSpecies]);
%             return
%         end
%     end
%     
%     % Temporary B particle matrix
%     bTemp = squeeze( particleMatrix(jSpecies, 1, veloRangeFull, radiusRangeFull) );
%     
%     % If the dimensions don't agree
%     if [numel(veloRangeFull) numel(radiusRangeFull)] ~= size(bTemp)
%         % Transpose array
%         bTemp = bTemp';
%     end
%     
%     % If there are non-collided B particles left
%     if any( bTemp > nMin, 'all' )
%         
%         % If there are enough non-collided B particles
%         if any( bTemp > C_AB, 'all' )
%             
%             % Remove non-collided B particles
%             bTemp(bTemp > C_AB) = bTemp(bTemp > C_AB) - C_AB(bTemp > C_AB);
%             particleMatrix(jSpecies, 1, veloRangeFull, radiusRangeFull) = ...
%                 reshape( bTemp, 1, 1, numel(veloRangeFull), numel(radiusRangeFull) );
%             
%             if any(particleMatrix < -nMin, 'all')
%                 disp('There are non-collided B particles, and there are more then collisions.');
%                 disp([iTime iRadius iVeloZero+iVelo iSpecies jSpecies]);
%                 return
%             end
%             
%         % If there are enough non-collided A particles
%         else
%             
%             % Remove remainder from collided B particles
%             particleMatrix(jSpecies, 2, veloRangeFull, radiusRangeFull) = ...
%                 particleMatrix(jSpecies, 2, veloRangeFull, radiusRangeFull) ...
%                 - reshape( C_AB - bTemp, 1, 1, numel(veloRangeFull), numel(radiusRangeFull) );
%             
%             % Remove all non-collided B particles
%             particleMatrix(jSpecies, 1, veloRangeFull, radiusRangeFull) = 0;
%             
%             if any(particleMatrix < -nMin, 'all')
%                 disp('There are non-collided B particles, but there are more collisions.');
%                 disp([iTime iRadius iVeloZero+iVelo iSpecies jSpecies]);
%                 return
%             end
%             
%         end
%     % If there are no non-collided B particles left
%     else
%         % Remove all collisions from collided B particles
%         particleMatrix(jSpecies, 2, veloRangeFull, radiusRangeFull) = ...
%             particleMatrix(jSpecies, 2, veloRangeFull, radiusRangeFull) ...
%             - reshape( C_AB, 1, 1, numel(veloRangeFull), numel(radiusRangeFull) );
%         
%         if any(particleMatrix < -nMin, 'all')
%             disp('There are no non-collided B particles.');
%             disp([iTime iRadius iVeloZero+iVelo iSpecies jSpecies]);
%             return
%         end
%     end
%     
%     % Scale the actually traveled distance by the difference in velocity
%     %   before and after collision
%     nRadiusTraveledA = round( 0.5.*(iVelo + ( velo(v_A) ./ velo(iVeloZero + iVelo) ) .* iVelo) );
%     
%     % Add collided A particles to collision matrix
%     collisionMatrix(iSpecies, 1, v_A, iRadius+iVelo) = ...
%         collisionMatrix(iSpecies, 1, v_A, iRadius+iVelo) ...
%         + reshape( sum(C_AB, 2), 1, 1, numel(veloRangeFull), 1);
%     
%     % Add collided B particles to collision matrix
%     collisionMatrix(jSpecies, 1, v_B, radiusRangeFull) = ...
%         collisionMatrix(jSpecies, 1, v_B, radiusRangeFull) ...
%         + reshape( C_AB, 1, 1, numel(veloRangeFull), numel(radiusRangeFull) );
%     
%     if any(collisionMatrix < -nMin, 'all')
%         disp('Negative collision matrix.');
%         disp([iTime iRadius iVeloZero+iVelo iSpecies jSpecies]);
%         return
%     end
% 
%     % Increment number of loops
%     NLoops = NLoops + 1;
%     
% end % For jSpecies

end % For iSpecies

end % For iVelo

end % For dir

end % For iRadius

%-------------------------------------------------------------------------------
% Update non-collided particles
%-------------------------------------------------------------------------------

% Subtract collided particles
particleMatrix = particleMatrix - subtractMatrix;

% Propagate non-collided particles
particleMatrix = updateMatrix( particleMatrix, nVelo, iVeloZero, keepParticleBool );

% Add collided particles back to particle matrix
particleMatrix(:, 2, :, :) = particleMatrix(:, 2, :, :) + collisionMatrix;

%-------------------------------------------------------------------------------
% Plot results
%-------------------------------------------------------------------------------

% Only for specific times
if (iAngle == 1) && any(iTime == plotTimes)
    
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
            
            % If non-collided and plasma species and velocity smoothing is off
            if ~veloSmoothNonColBool && (k == 1) && (iSpecies ~= 1)
                    % Total number of particles for this number of collisions
                    nParticlePerK = sum( particleMatrix(iSpecies, k, :, :), 'all' );

                    % Fit a normalized Gaussian curve to the number of particles
                    nParticleRadius = fitGaussian( radius, nParticleRadius, nMin );

                    % Multiply by number of particles
                    nParticleRadius = nParticleRadius .* nParticlePerK;
                else
            end
            
            % If plot smoothing is enabled and the density of the
            %   background gas is not plotted
            if smoothPlotBool && ~(plotDensityBool && (iSpecies == 1))
                % Smooth data
                nParticleRadius = ...
                    smoothdata( nParticleRadius, ...            % Data
                                2, ...                          % Dimension
                                'gaussian', ...                 % Type
                                round(nRadius * radiusMax) );   % Smoothing width
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
            if (iSpecies == 1) && plotDensityBool
                ylimMaxBG = max(nParticleRadiusSum ./ binVolume);
            else
                ylimMax = max(nParticleRadiusSum);
            end
        end
        
        % Set y-limit of plot
        if iSpecies == 1
            if plotDensityBool
                ylim([nMinBG ylimMaxBG]);
            else
                ylim([nMinBG ylimMax]);
            end
        else
            ylim([nMin ylimMax]);
        end
        
        % Set x-limit of plots
        xlim([0 max(radius)]);
        
        % Enable legend for last subplot
        if plotTimeIndex == numel(plotTimes)
            legend;
            legend('boxoff');
            legend('Location', 'northeast');
        end

        hold off;
        
    end % Species loop
end % Time if

% [DEBUG]
if debugBool
    % Reduced particle 2D matrices per particle type
    %   * Helpful for troubleshooting and to gain insight in particle propagation 
    bgMatrix = squeeze(sum(particleMatrix(1, :, :, :), 2));
    plasma1Matrix = squeeze(sum(particleMatrix(2, :, :, :), 2));
    plasma2Matrix = squeeze(sum(particleMatrix(3, :, :, :), 2));
    plasma3Matrix = squeeze(sum(particleMatrix(4, :, :, :), 2));
    
%     % Check for negative number of particles
%     if any( particleMatrix < -min([nMin nMinBG]), 'all')
%         error('Negative number of particles detected.');
%     end
    
    % If conservation of particle is enabled
    if keepParticleBool
        % Check if the total number of particles is conserved.
        if abs( sum( particleMatrix, 'all' ) - nParticleTotal ) > max([nMin nMinBG])
            error('The total number of particles is not conserved.');
        end
    end
end

% Show time
disp(['Simulation time: ' num2str(time(iTime)) ' s (Execution time: ' num2str(toc(durationTotalTimer)) ' s)']);

% Save execution times in array
executionTime(iTime) = toc(durationTotalTimer);

end % Time loop

end % Angle loop

% Calculate total number of collided and non-collided particles
nNonColTotal = sum( particleMatrix(:, 1, :, :), 'all' );
nColTotal = sum( particleMatrix(:, 2, :, :), 'all' );

% Calculate the final plasma and background gas temperatures
[plasmaTempEnd, bgTempEnd] = averageTemperature( particleMatrix, velo, mass );

% Get total excecution time [s]
durationTotal = toc(durationTotalTimer);

%% Save all open figures

if saveFiguresBool
    % Get all open figures
    figArray = findall(0,'Type','figure');
    
    % Loop through figures
    for iFig = 1 : numel(figArray)
        % Save figure in results folder with figure name followed by date / time
        savefig( figArray(iFig), ...
                 [ currentPath '\' folder '\' ...
                   get(figArray(iFig), 'Name') ' ' ...
                   figureLabel ' ', ...
                   datestr(now, 'yyyy-mm-dd HHMM') ...
                   saveFormat ] );
    end
end

%% Create config file

if createConfigBool
    createConfigFile( currentPath, folder, fileName, commentString, ... % File settings
        plotTimes, plotVeloDistBool, smoothPlotBool, plotDensityBool, plotLogBool, plot2DBool, saveFiguresBool, saveFormat, ... % Preferences
        collisionBool, nMin, nMinBG, kMax, negativeVeloBool, veloSmoothNonColBool, veloSmoothColBool, smoothWidth, keepParticleBool, ... % Computational restrictions
        time, angle, radius, velo, ... % Dimensional limits
        uc, ucRatio, atomUC, nAtomUC, nOxidePerElement, targetDensity, densityRatio, absorbedRatio, heatTarget, energyExcitation, ... % Material parameters
        bg, bgPressure, bgTemperature, bgDensity, ... % Background gas parameters
        spotWidth, spotHeight, ablationDepth, ablationVolume, laserFluence, energyLaser, ... % Laser parameters
        cosPowerFit, distributionType, initVeloDistWidth, ... % Model parameters
        nAtomAblated, nPlasmaTotal, nBGTotal, nNonColTotal, nColTotal, ... % Number of particles
        plasmaTempStart, plasmaTempEnd, bgTempStart, bgTempEnd, ... % Temperatures
        debugBool, durationTotal, nLoops ); % Technical                         
end
