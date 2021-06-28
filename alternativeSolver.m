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
plotVeloDistBool = false;

% Smooth plot results
smoothPlotBool = true;

% Plot density plot for background gas particles instead of number of particles
plotDensityBool = false;

% Plot figures in log scale
plotLogBool = false;

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
nMin    = 1E5;

% Maximum number of collisions
%   * Currently only non-collided (k=0) and collided (k=1) results get plotted
%   * Separating the plots into more different number of collisions is
%       supported, but not implemented at this moment.
kMax = 0;

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
radiusDelta = 5E-5;     % Radial step size          ( Default: 4-6E-5 )

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
if nMin > (nMin + 1)
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
for iSpecies = 1 : numel(uc)
    nAtomAblated = nAtomAblated + nUCAblated(iSpecies) .* sum([uc(iSpecies).AMOUNT]);
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

% Number of loops counter
NLoops = 0;

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
iVeloMaxAlt = round( ( veloInit(1).*(min(mass) - max(mass)) ...
                       + 2*max(mass)*veloInit(iVeloMax) ) ...
                     ./ ( min(mass) + max(mass) ) ./ veloDelta );

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
collisionMatrix = zeros(nSpecies, kMax+1, nVelo, nRadius);

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

% Recalculate number of particles each time step
nBGTotalNew = sum( particleMatrix(1, :, :, :), 'all' );
nPlasmaTotalNew = sum( particleMatrix(2:end, :, :, :), 'all' );

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
particleMatrix(particleMatrix < nMin) = 0;

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
        nLoops, nMin, nMin, iRadiusRange, radiusDelta, nRadius, velo, veloDelta, nVelo, ...
        iVeloRange, iVeloZero, binVolume, colCS, mass, nSpecies );
end

%-------------------------------------------------------------------------------
% New stuff
%-------------------------------------------------------------------------------

% Initialize number of collisions matrix
nCol = zeros(nSpecies, 1, nVelo, nRadius);

% Radii of particles
particleRadius = [sum([bg.RADIUS]) atomUC.RADIUS];

% Total number of collisions this time step
nColTotal = 0;

% For each possible particle speed
for iSpeed = iVeloMax-1 : -1 : 1
    
    % Reset number of collisions matrix
    nCol = nCol .* 0;
    
    % For each direction
    for iDir = [1 -1]

        % Velocity index
        iVelo = iSpeed * iDir;
        
        % For each plasma species
        for iSpecies = nSpecies : -1 : 2
            
            % For each possible position 
            for iRadius = 1 : nRadius
                
                % Skip radial bin if the number of A particles is below threshold
                if particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius) < nMin
                    continue
                end
                
                % If the particle is moving forwards (towards substrate)
                if iDir == 1
                    % All velocities a particle at the current velocity can collide with
                    iVeloCol = 1:iVeloZero+iVelo-1;

                    % The radial positions a particle with current velocity can collide
                    %   with within one time step
                    if iRadius+iVelo-1 > nRadius
                        iRadiusCol = iRadius:nRadius;
                    else
                        iRadiusCol = iRadius:iRadius+iVelo-1;
                    end
                % If the particle if moving backwards (towards target)
                elseif iDir == -1
                    % All velocities a particle at the current velocity can collide with
                    iVeloCol = iVeloZero+iVelo+1:nVelo;

                    % The radial positions a particle with current velocity can collide
                    %   with within one time step
                    if iRadius+iVelo+1 < 1
                        iRadiusCol = 1:iRadius;
                    else
                        iRadiusCol = iRadius+iVelo+1:iRadius;
                    end
                end
                
                % Skip radial bin if all B particles it is possible to
                %   collide with are below threahold
                if all( particleMatrix(1:iSpecies-1, 1, iVeloCol, iRadiusCol) < nMin, 'all' )
                    continue
                end
                
                % Number of A particles
                nA = particleMatrix(iSpecies, 1, iVeloZero+iVelo, iRadius);
                
                % Number of B particles
                nB = particleMatrix(1:iSpecies-1, 1, iVeloCol, iRadiusCol);
                nB = reshape( nB, numel(1:iSpecies-1), 1, numel(iVeloCol), numel(iRadiusCol) );
                
                % Relative velocities between A and B particles
                veloRelative = abs( velo(iVeloZero+iVelo) - velo(iVeloCol) );
                veloRelative = reshape( veloRelative, 1, 1, numel(iVeloCol) ) ;
                
                % Collisions cross-section of particles A and B
                sigma_AB = pi .* (particleRadius(iSpecies) + particleRadius(1:iSpecies-1)).^2;
                sigma_AB = reshape( sigma_AB, numel(1:iSpecies-1), 1 );
                
                % Initialize temporary number of collisions array
                nColTemp = zeros( numel(1:iSpecies-1), 1, numel(iVeloCol), numel(iRadiusCol) );
                
                % Loop through the traveled radial bins
                for jRadius = 1 : numel(iRadiusCol)
                    % Calculate the number of collisions within one timestep
                    nColTemp(:, :, :, jRadius) = nA .* nB(:, :, :, jRadius) .* sigma_AB .* timeDelta ...
                        .* veloRelative ./ binVolume(iRadiusCol(jRadius));
                    
                    % Update the remaining number of A particles
                    nA = nA - sum( nColTemp(:, 1, :, jRadius), 'all' );
                    
                    % Break loop if no particles remain
                    if nA < nMin
                        break
                    end
                    
                    % Increment number of loops
                    NLoops = NLoops + 1;
                end
                
                % Reshape to ensure the correct dimensions
                nColTemp = reshape( nColTemp, numel(1:iSpecies-1), 1, numel(iVeloCol), numel(iRadiusCol) );
                       
                % Remove collisions below threshold
                nColTemp(nColTemp < nMin) = 0;
                
                % Skip if no collisions remain
                if all( nColTemp == 0, 'all' )
                    continue
                end
                       
                % Sum of all collisions at this radius and velocity
                nColSum = sum( nColTemp, 'all' );
                
                % Add A particle collisions to number of collision matrix
                nCol(iSpecies, 1, iVeloZero+iVelo, iRadius) = ...
                    nCol(iSpecies, 1, iVeloZero+iVelo, iRadius) ...
                    + 0.5*nColSum;
                
                % Add B particle collisions to number of collisions matrix
                nCol(1:iSpecies-1, 1, iVeloCol, iRadiusCol) = ...
                     nCol(1:iSpecies-1, 1, iVeloCol, iRadiusCol) ...
                     + 0.5.*nColTemp;
                 
                % Velocity of particle A after collision
                veloA = ( ( mass(iSpecies) - mass(1:iSpecies-1) ) * velo(iVeloZero+iVelo) ...
                        + 2 * mass(1:iSpecies-1) .* velo(iVeloCol) ) ...
                      ./ ( mass(iSpecies) + mass(1:iSpecies-1) );

                % Velocity of particles B after collision
                veloB = ( ( mass(1:iSpecies-1) - mass(iSpecies) ) .* velo(iVeloCol) ...
                        + 2 * mass(iSpecies) .* velo(iVeloZero+iVelo) ) ...
                      ./ ( mass(iSpecies) + mass(1:iSpecies-1) );

                % Velocity indices after collision
                veloA = iVeloZero + round( veloA ./ veloDelta );
                veloB = iVeloZero + round( veloB ./ veloDelta );

                % Prevent index out-of-bounds error
                %   * The code should be written such that this error is impossible to occur
                %   * So room for improvement here
                veloA(veloA > nVelo) = nVelo;
                veloA(veloA < 1) = 1;
                veloB(veloB > nVelo) = nVelo;
                veloB(veloB < 1) = 1;

                % Reshape to the right dimensions
                veloA = reshape( veloA, numel(1:iSpecies-1), numel(iVeloCol) );
                veloB = reshape( veloB, numel(1:iSpecies-1), numel(iVeloCol) );
                
                % Loop through other species
                for jSpecies = 1 : iSpecies-1
                    % Add A particle collisions to number of collision matrix
                    collisionMatrix(iSpecies, 1, veloA(jSpecies, :), iRadiusCol(jRadius)) = ...
                        collisionMatrix(iSpecies, 1, veloA(jSpecies, :), iRadiusCol(jRadius)) ...
                        + reshape( sum( nColTemp(jSpecies, 1, :, :), 4 ), 1, 1, numel(iVeloCol), 1);

                    % Add B particle collisions to number of collisions matrix
                    collisionMatrix(jSpecies, 1, veloB(jSpecies, :), iRadiusCol) = ...
                         collisionMatrix(jSpecies, 1, veloB(jSpecies, :), iRadiusCol) ...
                         + reshape( nColTemp(jSpecies, 1, :, :), 1, 1, numel(iVeloCol), numel(iRadiusCol) );
                
                    % Increment number of loops
                    NLoops = NLoops + 1;
                end % jSpecies
                
            end % iRadius
        end % iSpecies
    end % iDir
    
    % Limit number of collisions by number of particles
    nCol(nCol > particleMatrix) = particleMatrix(nCol > particleMatrix);
    
    % Add number of collision for this velocity to the total number of
    %   collisions this thime step
    nColTotal = nColTotal + sum( nCol, 'all' );
    
%     if abs( sum( collisionMatrix, 'all' ) - nColTotal ) > nMin
%         disp(sum( collisionMatrix, 'all' ));
%         disp(nColTotal);
%         return
%     end
    
    % Subtract collided perticles
    particleMatrix = particleMatrix - nCol;
    
end % iSpeed

collisionMatrix = ( collisionMatrix ./ sum( collisionMatrix, 'all' ) ) .* nColTotal;

%-------------------------------------------------------------------------------
% Update non-collided particles
%-------------------------------------------------------------------------------

% % Subtract collided particles
% particleMatrix = particleMatrix - subtractMatrix;

% Propagate non-collided particles
particleMatrix = updateMatrix( particleMatrix, nVelo, iVeloZero, keepParticleBool );

particleMatrix = particleMatrix + collisionMatrix;

% Add collided particles back to particle matrix
% particleMatrix(:, 2, :, :) = particleMatrix(:, 2, :, :) + collisionMatrix;

%-------------------------------------------------------------------------------
% Plot results
%-------------------------------------------------------------------------------

% Only for specific times
if (iAngle == 1)
    plotParticles( particleMatrix, plotTimes, iTime, time, atomUC, ...
        nMin, radius, binVolume, smoothPlotBool, plotDensityBool, plotLogBool );
end

% [DEBUG]
if debugBool
    % Reduced particle 2D matrices per particle type
    %   * Helpful for troubleshooting and to gain insight in particle propagation 
    bgMatrix = squeeze(sum(particleMatrix(1, :, :, :), 2));
    plasma1Matrix = squeeze(sum(particleMatrix(2, :, :, :), 2));
    plasma2Matrix = squeeze(sum(particleMatrix(3, :, :, :), 2));
    plasma3Matrix = squeeze(sum(particleMatrix(4, :, :, :), 2));
    
    % Check for negative number of particles
    if any( particleMatrix < -min([nMin nMin]), 'all')
        error('Negative number of particles detected.');
    end
    
%     % If conservation of particle is enabled
%     if keepParticleBool
%         % Check if the total number of particles is conserved.
%         if abs( sum( particleMatrix, 'all' ) - nParticleTotal ) > max([nMin nMin])
%             error('The total number of particles is not conserved.');
%         end
%     end
end

% Show time
disp(['Simulation time: ' num2str(time(iTime)) ' s (Execution time: ' num2str(toc(durationTotalTimer)) ' s)']);

% Save execution times in array
executionTime(iTime) = toc(durationTotalTimer);

end % Time loop

end % Angle loop

% Calculate total number of collided and non-collided particles
nNonColTotal = sum( particleMatrix(:, 1, :, :), 'all' );
% nColTotal = sum( particleMatrix(:, 2, :, :), 'all' );
nColTotal = 0;

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
        collisionBool, nMin, nMin, kMax, negativeVeloBool, veloSmoothNonColBool, veloSmoothColBool, smoothWidth, keepParticleBool, ... % Computational restrictions
        time, angle, radius, velo, ... % Dimensional limits
        uc, ucRatio, atomUC, nAtomUC, nOxidePerElement, targetDensity, densityRatio, absorbedRatio, heatTarget, energyExcitation, ... % Material parameters
        bg, bgPressure, bgTemperature, bgDensity, ... % Background gas parameters
        spotWidth, spotHeight, ablationDepth, ablationVolume, laserFluence, energyLaser, ... % Laser parameters
        cosPowerFit, distributionType, initVeloDistWidth, ... % Model parameters
        nAtomAblated, nPlasmaTotal, nBGTotal, nNonColTotal, nColTotal, ... % Number of particles
        plasmaTempStart, plasmaTempEnd, bgTempStart, bgTempEnd, ... % Temperatures
        debugBool, durationTotal, nLoops ); % Technical                         
end
