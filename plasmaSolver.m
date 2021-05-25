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
debugBool = false;

% Enable / disable timers to measure code performance
timerBool = true;

% Create a configuration file containing simulation settings
createConfigBool = true;

% Plot the initial velocity distribution
plotInitVeloDistBool = true;

% Save all figures
saveFiguresBool = false;

% Plot all angles (true) or just the centre of the plasma (false)
plot2DBool = false;

%-------------------------------------------------------------------------------
% Computational restrictions
%-------------------------------------------------------------------------------

% Minimal number of particles per bin
%   * Best way to increase performance, increasing the value to much will
%     result in non-sensicle results.
nMin = 1E8;

% Maximum number of collisions
%   * Currently only non-collided (k=0) and collided (k=1) results get plotted
kMax = 1;

% Switch velocity smoothing for non-collided particles on or off (true / false)
%   * Signifiacnt impact on performance
%   * Prevents quantization error
%   * Conserves number of particles
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
timeMax     = 10E-6;    % Simulation run time ( Default: 5-20E-6 )
timeDelta   = 1E-7;     % Time step duration  ( Default: 1E-7    )

% Angular limits [deg]
angleMax    = 90;       % End angle         ( Default: 90 )
angleDelta  = 3;        % Angular step size ( Default: 3  )

% Radial limits [m]
radiusMax   = 0.05;     % Target-substrate distance ( Default: 0.05 )
radiusDelta = 2E-5;     % Radial step size          ( Default: 6E-5 ) 

% Maximum velocity [m / s]
%   * Increase if target contains atoms lighter than oxygen (e.g. lithium)
%   * Default: 3E4
veloMax = 3E4;

%-------------------------------------------------------------------------------
% Material parameters
%-------------------------------------------------------------------------------

% Unit cell
%   * Preface with UC.
%   * Any mixture of crystals is possible, eg. [UC.TiO2 UC.Li2O]
%   * To add new materials or to view material propertie, see UnitCells.m
%   Options: TiO2, SrTiO3, Li2O, LiO2, LiTi2O4, Li2TiO3, Li4Ti5O12
uc = [UC.TiO2 UC.Li2O];

% Mixture ratio of species in target
%   * Ignore if target contains only one component
%   * Should be of same length as unit cell array
ucRatio = [1 2];

% Target density [kg / m^3]
%   * If set to 0, a perfect single crystal target is assumed
%   * For most accurate results measure per target
%   * Default: 0
targetDensity = 0;

% Ratio of laser energy that is absorbed by the target
%   * Determine per target for most accurate results
%   * Default: 0.6
absorbedRatio = 0.6;

% The number of oxides that can form per element
%   * Example: From a TiO2 target TiO and TiO2 can form during propagation,
%       so nOxidesPerElement is 2.
%   * Default: 2
nOxidePerElement = 2;

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
bgPressure      = 0.1E2; % 0.02E2;

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
spotHeight      = 2.1E-3;

% Ablation depth into the target for a single pulse [m]
%   * For most accurate results calculate per target
%   * Default: 100E-9
ablationDepth   = 100E-9;

% Laser fluence [J / cm^2]
%   * Default: 2.0
laserFluence    = 2.0;

%-------------------------------------------------------------------------------
% Model parameters
%-------------------------------------------------------------------------------

% Fit parameter of the cosines power in the angular plasma particle distibution
%   function
%  * Default: 25
cosPowerFit = 25;

% Initial velocity distribution width
%   * Default: 1500
initVeloDistWidth = 1500;

%% Calculations
% Everything below is calculated automatically based on user input above.

%-------------------------------------------------------------------------------
% Timer
%-------------------------------------------------------------------------------

if timerBool
    % Time the entire script
    durationTotalTimer = tic;
end

% Initialize measured times
%   * To prevent createConfigFile throwing an error.
durationTotal = 0;

%-------------------------------------------------------------------------------
% Material calculations
%-------------------------------------------------------------------------------

% % Number of elements in target
% nElements = numel(uc.ELEMENTS);

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
    ucRatio = ucRatio / sum(ucRatio);
    
    % Weighted average single crystal density
    scDensity = sum(([uc.DENSITY] .* ucRatio) ./ sum(ucRatio));
    
    % Weighted average unit cell volume
    ucVolume = sum(([uc.VOLUME] .* ucRatio) ./ sum(ucRatio));
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

% Number of possible species
nSpecies = calculateNSpecies(atomUC, bg, nOxidePerElement);

% Get particle masses
mass = [sum([bg.MASS]) atomUC.MASS];

% Get metals present in target
metalUC = atomUC([atomUC.NUMBER] ~= 8);

% If a non-oxygen element is present in the target and there is oxygen
%   present in the background gas
if (numel(metalUC) > 1) && sum([bg.NUMBER] == 8)
    % Add oxides to mass array
    %   * Assumes oxides of form MO, MO2, MO3, etc.
    %   * This is only correct for some material, if a material forms
    %       different types of oxides it has to be added manually (e.g. Li2O)
    for iOxide = 1 : nOxidePerElement
        mass = [ mass ([metalUC.MASS] + iOxide*PT.O.MASS) ];
    end
end

% Get radii of plasma species
colCS = [atomUC.RADIUS];

if (numel(metalUC) > 1) && sum([bg.NUMBER] == 8)
    % Add radii of oxides
    %   * Assumes oxides of form MO, MO2, MO3, etc.
    %   * This is only correct for some material, if a material forms
    %       different types of oxides it has to be added manually (e.g. Li2O)
    for iOxide = 1 : nOxidePerElement
        colCS = [ colCS ([metalUC.RADIUS] + iOxide*PT.O.RADIUS) ];
    end
end

% Calculate collisions cross-sections
colCS = pi .* (sum([bg.RADIUS]) + colCS) .^2;

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
nUCAblated = ablationVolume ./ [uc.VOLUME];

% Number of ablated atoms
nAtomAblated = (ablationVolume / ucVolume) .* nAtomUC;

%% Axis creation

% Temporal axis
time    = 0 : timeDelta : timeMax - timeDelta;
nTime   = numel(time);

% Angular axis
angle   = 0 : angleDelta : angleMax - angleDelta;

% Set number of angles
if plot2DBool
    nAngle  = numel(angle);
else
    nAngle = 1;
end

% Radial axis
radius  = 0 : radiusDelta : radiusMax - radiusDelta;
nRadius = numel(radius);

% Velocity array
veloDelta   = radiusDelta / timeDelta;
velo        = 0 : veloDelta : veloMax - veloDelta;
nVelo       = numel(velo);

%% Pre-allocations

% Particle matrix 
particleMatrix = zeros(nSpecies, kMax+1, nVelo, nRadius);

% Matrix holding new positions of collided particles
%   * Reset every time step
collisionMatrix = zeros(nSpecies, 1, nVelo, nRadius);

% Matrix holding all collision rates
%   * Calculated every time step
%   * Prevents having to calculate the collision rate millions of times when
%       when called, instead calculate the collision rate for every position at
%       once each time step
collisionRate = zeros(nSpecies-1, 1, nVelo, nRadius);

% Matrices holding indices of new velocity of particles after collision
iNewVeloMatrix = zeros(nVelo, nSpecies-1, nVelo);   % Plasma particles
iNewVeloBGMatrix = zeros(nVelo, nSpecies-1, nVelo); % Background gas particles

%% Pre-calculations

% Matrix containing all possible velocity weights for collision calculation
%   * Prevents having to calculate many times, instead just look up value
%       from matrix
veloWeight = (velo' - velo) ./ (velo' + velo);

% Loop through plasma species
for iSpecies = 2 : nSpecies
    % Calculate new velocity of plasma species after collision with a
    %   background gas particle
    iNewVeloMatrix(:, iSpecies, :) = round( (velo'.*(mass(iSpecies) - mass(1)) ...
        + 2*mass(1).*velo) ./ (mass(iSpecies) + mass(1)) ./ veloDelta );
    
    % Calculate new velocity of background gas particle after collision with
    %   plasma species
    iNewVeloBGMatrix(:, iSpecies, :) = round( (velo'.*(mass(1) - mass(iSpecies)) ...
        + 2*mass(iSpecies).*velo) ./ (mass(iSpecies) + mass(1)) ./ veloDelta );
end

% Prevent index out-of-bounds error
%   * Particles are not allowed to move backwards, so every velocity below zero
%       is set to zero
%   * Particles cannot exceed maximum velocity
iNewVeloMatrix(iNewVeloMatrix < 1) = 1;
iNewVeloBGMatrix(iNewVeloBGMatrix < 1) = 1;
iNewVeloMatrix(iNewVeloMatrix > nVelo) = nVelo;
iNewVeloBGMatrix(iNewVeloBGMatrix > nVelo) = nVelo;

% Angular plasma particle distribution
nParticleAngle = angularDistribution( angle,        ...
                                      radiusDelta,  ...
                                      cosPowerFit,  ...
                                      sum(nAtomAblated) );

%% Main program

for iAngle = 1 : nAngle
%% Calculations per angle

% Skip angle if there are not enough particles present
if nParticleAngle(iAngle) < nMin
    continue
end

% Reset particle matrix
if iAngle > 1
    particleMatrix = particleMatrix.*0;
end

% Pre-allocate background gas particle matrix and fill matrix based on
%   calculated density, also retrieve bin volumes per radial bin
[particleMatrix(1, 1, 1, :), binVolume] = fillBGMatrix( bgDensity,    ...
                                                        radius,       ...
                                                        angle,        ...
                                                        iAngle,       ...
                                                        radiusDelta );
                                        
% Total number of background particles at this angle
nBGTotal = sum( particleMatrix(1, 1, 1, :) );

% Total number of plasma particles at this angle
nPlasmaTotal = nParticleAngle(iAngle);

% Calculate the initial particle velocity distribution
nParticleVeloInit = initialVelocityDistribution( uc,                ...
                                                 ucRatio,           ...
                                                 atomUC,            ...
                                                 ablationVolume,    ...
                                                 velo,              ...
                                                 initVeloDistWidth, ...
                                                 absorbedRatio,     ...
                                                 energyLaser,       ...
                                                 heatTarget,        ...
                                                 energyExcitation,  ...
                                                 plotInitVeloDistBool );

return                                           

% Fill into the plasma matrix
particleMatrix(2:3, 1, :, 1) = reshape(nParticleVeloInit(:, 1:2)', nElements, 1, nVelo, 1);

nLoops = 0;
errorCount = 0;

nMinSum = nMin; % nMin * (nSpecies - 1) * (kMax+1);

colorArray = ['b' 'r' 'k'];
plotArray = {'Uncollided', 'Collided', 'Total'};

plotIndex = 0;

for iTime = 1 : 1 % nTime
%% Calculations per time step

% Reset collision matrix
collisionMatrix = collisionMatrix.*0;

% Total number of particles
particleTotal = sum( particleMatrix, 'all' );

% Skip first time step for smoothing
if iTime > 1
    if veloSmoothNonColBool
        % Smooth velocities of non-collided particles
        particleMatrix(:, 1, 2:end, :) = smoothdata( particleMatrix(:, 1, 2:end, :), 3, 'gaussian', smoothWidth );
    end
    
    if veloSmoothColBool
        % Smooth velocities of collided particles
        particleMatrix(:, 2:end, 2:end, :) = smoothdata( particleMatrix(:, 2:end, 2:end, :), 3, 'gaussian', smoothWidth );
    end
end

% Remove all particles below threshold
particleMatrix(particleMatrix < nMin) = 0;

% Normalize number of particles and add removed particles back to filled bins
particleMatrix = ( particleMatrix ./ sum(particleMatrix, 'all') ) .* particleTotal;

% Calculate collision rate matrix
collisionRate = ( sum(particleMatrix(1, :, :, :), 2) ./ reshape(binVolume, 1, 1, 1, nRadius) ) ...
                .* radiusDelta .* colCS;

timerA = tic;
for iRadius = (nRadius - 1) : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

for iVelo = nVelo : -1 : 2
%% Calculations per non-oxidizing plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Only loop through velocities that traverse at least one radial bin per
%   time step

% If number of plasma particles is zero
if sum(particleMatrix(2:end, :, iVelo, iRadius), 2) == 0
    continue
end

% Remaining number of particles
nPlasmaTemp = particleMatrix(2:end, :, iVelo, iRadius);

% Restrict number of traveled bins to the total number of
%   radial bins
if ( iRadius + iVelo - 1 ) > nRadius
    nRadiusDelta = nRadius - iRadius;
else
    nRadiusDelta = iVelo - 1;
end

for jRadius = 0 : (nRadiusDelta - 1)
%% Calculations per traversed radial bin
% Loop from current bin to one bin before final bin
% Collisions occur first in the bin closest to starting bin

% Break loop if no more particles are left to collide
if sum(nPlasmaTemp, 'all') < nMinSum
    break
end

% Index of the current radius: starting bin + traveled bin
thisRadius = iRadius + jRadius;

for iVeloBG = 1 : (iVelo - 1)
%% Calculations per background velocity bin
% Only include background velocities smaller than the plasma velocity
% From slowest to highest as collisions with slow moving background
%   particles are more probable

timerC = tic;
% Break loop if no more particles are left to collide
if sum(nPlasmaTemp, 'all') < nMinSum
    break
end

% If number of background particles is below threshold
if sum(particleMatrix(1, :, iVeloBG, thisRadius), 2) == 0
    continue
end

% Number of collisions array
nCol = veloWeight(iVelo, iVeloBG) .* collisionRate(:, 1, iVeloBG, thisRadius) .* nPlasmaTemp;

% Skip iteration if the number of collisions is smaller than the threshold
if sum(nCol, 'all') < nMinSum
    continue
end

% If the number of collisions is greater than the number of background particles
if sum(nCol, 'all') > sum(particleMatrix(1, :, iVeloBG, thisRadius), 2)
    % Normalize number of collisions to number of background particles
    %   (minus 1 to avoid faulty comparisons)
    nCol = (nCol ./ sum(nCol, 'all')) .* (sum(particleMatrix(1, :, iVeloBG, thisRadius), 2) - 1);
end

% Update number of non-collided plasma particles
nPlasmaTemp = nPlasmaTemp - nCol;
timeC = toc(timerC);

%--------------------------------------------------------------------------
% Calculate new radial positions after collision
%--------------------------------------------------------------------------

timerD = tic;
% Reshape indice arrays
iNewVelo = iNewVeloMatrix(iVelo, :, iVeloBG);
iNewVeloBG = iNewVeloBGMatrix(iVelo, :, iVeloBG);

% New radial index of collided plasma particles
nNewRadiusDelta = round( velo(iNewVelo) .* (nRadiusDelta - jRadius) ...
                         .* timeDelta / (nRadiusDelta * radiusDelta) );
                       
% New radial index of collided background particles
nNewRadiusDeltaBG = round( velo(iNewVeloBG) .* (nRadiusDelta - jRadius) ...
                           .* timeDelta / (nRadiusDelta * radiusDelta) );

% Prevent index out-of-range error
nNewRadiusDelta(thisRadius + nNewRadiusDelta > nRadius) = ...
    nRadius - thisRadius;
nNewRadiusDeltaBG(thisRadius + nNewRadiusDeltaBG > nRadius) = ...
    nRadius - thisRadius;
timeD = toc(timerD);

%-------------------------------------------------------------------------------
% Update matrices
%-------------------------------------------------------------------------------

timerE = tic;
% Remove plasma particles from initial position before collision
particleMatrix(2:end, :, iVelo, iRadius) = particleMatrix(2:end, :, iVelo, iRadius) - nCol;

% Remove background particles from initial position before collision
particleMatrix(1, :, iVeloBG, thisRadius) = particleMatrix(1, :, iVeloBG, thisRadius) ...
    - ( particleMatrix(1, :, iVeloBG, thisRadius) ./ sum(particleMatrix(1, :, iVeloBG, thisRadius)) ) * sum(nCol, 'all');

nColSum = sum(nCol, 2);

for iSpecies = 1 : 4 % (nSpecies - 1)
    % Add plasma particles to new position after collision
    collisionMatrix(iSpecies+1, 1, iNewVelo(iSpecies), thisRadius + nNewRadiusDelta(iSpecies)) = ...
        collisionMatrix(iSpecies+1, 1, iNewVelo(iSpecies), thisRadius + nNewRadiusDelta(iSpecies)) + nColSum(iSpecies);

    % Add background particles to new position after collision
    collisionMatrix(1, 1, iNewVeloBG(iSpecies), thisRadius + nNewRadiusDeltaBG(iSpecies)) = ...
        collisionMatrix(1, 1, iNewVeloBG(iSpecies), thisRadius + nNewRadiusDeltaBG(iSpecies)) + nColSum(iSpecies);
end
timeE = toc(timerE);

nLoops = nLoops + 1;

end % Background velocity loop

end % Traveled distance loop

end % Velocity loop

end % Radius loop
timeA = toc(timerA);

timerB = tic;
particleMatrix = updateMatrix(particleMatrix, nVelo, nSpecies, kMax+1, ...
                                keepParticleBool );

particleMatrix(:, 2, :, :) = particleMatrix(:, 2, :, :) + collisionMatrix;
timeB = toc(timerB);

% Only for specific times
if (iTime ==  6) || (iTime == 11) || (iTime == 21) || (iTime == 31) || (iTime == 61)
    
    plotIndex = plotIndex + 1;
    
    % Loop through plasma species
    for iSpecies = 1 : 3        
        figure(iSpecies);
        subplot(1, 5, plotIndex);
        title([num2str(time(iTime), 3) ' s']);
        hold on;
        
        % Initialize sum of all collisions array
        nParticleRadiusSum = zeros(1, nRadius);
        
        for k = 1 : 2 % kMax+1
            % Calculate total number of particles per radial bin
            nParticleRadius = nParticlesPerRadius( radius, particleMatrix(iSpecies, k, :, :) );

            if k == 1 && iSpecies ~= 1
                % Total number of particles for this number of collisions
                nParticlePerK = sum( particleMatrix(iSpecies, k, :, :), 'all' );

                % Fit a normalized Gaussian curve to the number of particles
                nParticleRadius = fitGaussian( radius, nParticleRadius, nMin )';

                % Multiply by number of particles
                nParticleRadius = nParticleRadius .* nParticlePerK;
            else
                % Smooth data
                nParticleRadius = smoothdata( nParticleRadius, 2, 'gaussian', nRadius / 100 );
            end
            
            % Add to sum of all collisions array
            nParticleRadiusSum = nParticleRadiusSum + nParticleRadius;
            
            if iSpecies == 1
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
        end
       
        if iSpecies == 1
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
        
        xlim([0 0.05]);
        legend;
        hold off;
        
        veloMeanSquared = zeros(1, nRadius);
        tempMean = veloMeanSquared;
        
        figure(3+iSpecies);
        subplot(1, 2, 1);
        hold on;
        for iAtom = 1 : nRadius
            veloMeanSquared(iAtom) = sum(reshape(sum(particleMatrix(iSpecies, :, :, iAtom), 2), 1, nVelo) .* (velo.^2)) / sum(particleMatrix(iSpecies, :, :, iAtom), 'all');
            
            if iSpecies == 1
                tempMean(iAtom) = (massBG * veloMeanSquared(iAtom)) / (3 * C.BOLTZMANN);
            else
                tempMean(iAtom) = (mass(iSpecies) * veloMeanSquared(iAtom)) / (3 * C.BOLTZMANN);
            end
        end
        plot( radius, sqrt(veloMeanSquared), ...
              'DisplayName', [num2str(time(iTime), 3) ' s'] );
        hold off;
        subplot(1, 2, 2);
        hold on;
        plot( radius, tempMean, ...
              'DisplayName', [num2str(time(iTime), 3) ' s'] );
        hold off;
        
    end
end

end % Time loop

end % Angle loop

durationTotal = toc(durationTotalTimer);

%% Create config file

if createConfigBool
    createConfigFile(currentPath, folder, fileName, commentString, ...          % File settings
        saveFiguresBool, plot2DBool, ...                                        % Preferences
        nMin, kMax, veloSmoothNonColBool, veloSmoothColBool, smoothWidth, keepParticleBool, ... % Computational restrictions
        time, angle, radius, velo, ...                                          % Dimensional limits
        uc, ucRatio, targetDensity, absorption, ...                             % Material parameters
        bg, bgPressure, bgTemperature, bgDensity, ...                           % Background gas parameters
        spotWidth, spotHeight, ablationDepth, ablationVolume, laserFluence, energyLaser, heatTarget, ... % Laser parameters
        cosPowerFit, initVeloDistWidth, ...                                      % Model parameters
        nPlasmaTotal, nBGTotal, ...                                             % Number of particles
        durationTotal );                                                        % Timers                          
end
