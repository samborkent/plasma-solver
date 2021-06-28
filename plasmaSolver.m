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

% For best readability increasing the right-hand text limit and the comment
%   formatting limit column widths to 80 is adviced.

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
    'Code demonstration' ...
    ];

% Add label to saved figures
figureLabel = 'TiO2 SrTiO3 0.02mbar';

%-------------------------------------------------------------------------------
% Preferences
%-------------------------------------------------------------------------------

% Create a configuration file containing simulation settings
createConfigBool = true;

% Enable / disable debug mode
%   * Checks conservation of number of particles
%   * Checks for negative number of particles
debugBool = false;

%-------------------------------------------------------------------------------
% Plot settings
%-------------------------------------------------------------------------------

% Times to plot [us]
plotTimes = [0.5 1.5 3 6 9];

% Plot the initial velocity distribution
plotVeloDistBool = true;

% Plot the propagation of particles
plotResults = true;

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
saveFiguresBool = true;

% Format to save figures in
%   * Ex.: '.fig', '.svg', '.png', etc.
saveFormat = {'.fig' '.svg'};

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
collisionBool = true;

% Minimal number of particles per bin
%   * Best way to improve performance
%   * Increases quantization error
%   * Don't increase the values too much, otherwise too many particle bins will
%       be skipped.
%   * Experiment to maximize performance / quality of results
%   * Default: 1E5 / 1E6 (Limit seems to be ~1E8)
nMin    = 1E9;     % Plasma threshold

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
timeMax     = 7E-6;    % Simulation run time        ( Default: 5-20E-6 )
timeDelta   = 1E-7;     % Time step duration        ( Default: 1E-7    )

% Angular limits [deg]
angleMax    = 90;       % End angle                 ( Default: 90 )
angleDelta  = 3;        % Angular step size         ( Default: 3  )

% Radial limits [m]
%   * Decreasing radiusDelta increases the resolution of the model
%   * Has significant impact on performance
radiusMax   = 0.05;     % Target-substrate distance ( Default: 0.05 )
radiusDelta = 10E-5;     % Radial step size          ( Default: 4-6E-5 )

%-------------------------------------------------------------------------------
% Material parameters
%-------------------------------------------------------------------------------

% Unit cell
%   * Preface with UC.
%   * Any mixture of crystals is possible, eg. [UC.TiO2 UC.Li2O]
%   * To add new materials or to view material propertie, see UnitCells.m
%   Options: TiO2, SrTiO3, LiTi2O4, Li4Ti5O12, etc.
uc = [UC.TiO2 UC.SrTiO3];

% Mixture ratio of species in target
%   * Ignore if target contains only one component
%   * Should be of same length as unit cell array
ucRatio = [1 2];

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

% Ratio of laser energy (248 nm) that is absorbed by the target
%   * Determine per target for most accurate results
%   * TiO2: 0.66, SrTiO3: 0.77
%   * Default: 0.7
absorbedRatio = 0.7;

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
% Times
%-------------------------------------------------------------------------------

% Convert to indices
plotTimes = round( (plotTimes.*1E-6) ./ timeDelta );

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

% Calculate the collision cross-sections
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
%   * The 1E4 term convert cm^-2 to m^-2
energyLaser = (laserFluence * 1E4) * (spotWidth * spotHeight);

%-------------------------------------------------------------------------------
% Atom calculations
%-------------------------------------------------------------------------------

% Number of ablated unit cells
%   * Calculated as the ablation volume, divided by the weighted average of the
%       unit cell volumes present in the target, multiplied by the density
%       ratio.
nUCAblated = (ablationVolume / sum([uc.VOLUME] .* ucRatio)) ...
             .* ucRatio .* densityRatio;

% Number of ablated atoms
%   * Calculated as the number of unit cells for each target component,
%       multiplied by the number of atoms in that unit cell.
nAtomAblated = 0;
for iSpecies = 1 : numel(uc)
    nAtomAblated = nAtomAblated ...
        + nUCAblated(iSpecies) .* sum([uc(iSpecies).AMOUNT]);
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
nTime   = numel(time);  % Number of time steps in the model

% Array to store the execution time, can be plotted against the
executionTime = zeros(1, nTime);

% Radial axis
radius  = 0 : radiusDelta : radiusMax - radiusDelta;
nRadius = numel(radius);    % Number of spatial steps in the model
iRadiusRange = 1 : nRadius;

% Velocity step size
veloDelta = round( radiusDelta / timeDelta );

% Initial velocity axis
veloInit = round(0 : veloDelta : 20E4);

%-------------------------------------------------------------------------------
% Angular distribution
%-------------------------------------------------------------------------------

% Angular plasma particle distribution
nPlasmaAngle = angularDistribution( angle,        ...
                                    angleDelta,   ...
                                    radiusDelta,  ...
                                    cosPowerFit,  ...
                                    nAtomAblated );

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

% Remove all particles below threshold
particleMatrix(particleMatrix < nMin) = 0;

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
% Update non-collided particles
%-------------------------------------------------------------------------------

% Propagate non-collided particles
particleMatrix = updateMatrix( particleMatrix, nVelo, iVeloZero, keepParticleBool );

% Add collided particles back to particle matrix
particleMatrix(:, 2, :, :) = particleMatrix(:, 2, :, :) + collisionMatrix;

%-------------------------------------------------------------------------------
% Plot results
%-------------------------------------------------------------------------------

% Only for the first angle
if iAngle == 1 && plotResults
    plotParticles( particleMatrix, plotTimes, iTime, time, atomUC, ...
        nMin, radius, binVolume, smoothPlotBool, plotDensityBool, plotLogBool );
end

%-------------------------------------------------------------------------------
% DEBUG
%-------------------------------------------------------------------------------

if debugBool
    % Reduced particle 2D matrices per particle type
    %   * Helpful for troubleshooting and to gain insight in particle propagation 
    bgMatrix = squeeze(sum(particleMatrix(1, :, :, :), 2));
    plasma1Matrix = squeeze(sum(particleMatrix(2, :, :, :), 2));
    plasma2Matrix = squeeze(sum(particleMatrix(3, :, :, :), 2));
    plasma3Matrix = squeeze(sum(particleMatrix(4, :, :, :), 2));
    
    % Check for negative number of particles
    if any( particleMatrix < -nMin, 'all')
        error('Negative number of particles detected.');
    end
    
    % If conservation of particle is enabled
    if keepParticleBool
        % Check if the total number of particles is conserved.
        if abs( sum( particleMatrix, 'all' ) - nParticleTotal ) > nMin
            error('The total number of particles is not conserved.');
        end
    end
end

%-------------------------------------------------------------------------------
% Excecution time
%-------------------------------------------------------------------------------

% Show time
disp(['Simulation time: ' num2str(time(iTime)) ' s (Execution time: ' num2str(toc(durationTotalTimer)) ' s)']);

if debugBool
    % Save execution times in array
    executionTime(iTime) = toc(durationTotalTimer);
end

end % Time loop

end % Angle loop

%% End calculations

% Calculate total number of collided and non-collided particles
nNonColTotal = sum( particleMatrix(:, 1, :, :), 'all' );
nColTotal = sum( particleMatrix(:, 2, :, :), 'all' );

% Calculate the final plasma and background gas temperatures
[plasmaTempEnd, bgTempEnd] = averageTemperature( particleMatrix, velo, mass );

%-------------------------------------------------------------------------------
% Excecution time
%-------------------------------------------------------------------------------

% Get total excecution time [s]
durationTotal = toc(durationTotalTimer);

% Plot the execution time versus the simulation time to visualize code performance.
if debugBool
    figure(nSpecies+2);
    scatter(time, executionTime);
    xlabel('Simulation time [s]');
    ylabel('Execution time [s]');
end

%% Save all open figures

if saveFiguresBool
    % Get all open figures
    figArray = findall(0,'Type','figure');
    
    % Loop through figures
    for iFig = 1 : numel(figArray)
        for iFormat = 1 : numel(saveFormat)
            % Save figure in results folder with figure name followed by
            %   date and time
            if strcmp( char( saveFormat(iFormat) ), '.fig' )
                savefig( figArray(iFig), ...
                         [ currentPath '\' folder '\' ...
                           get(figArray(iFig), 'Name') ' ' ...
                           figureLabel ' ', ...
                           datestr(now, 'yyyy-mm-dd HHMM') ...
                           char( saveFormat(iFormat) ) ] );
            else
                saveas( figArray(iFig), ...
                         [ currentPath '\' folder '\' ...
                           get(figArray(iFig), 'Name') ' ' ...
                           figureLabel ' ', ...
                           datestr(now, 'yyyy-mm-dd HHMM') ...
                           char( saveFormat(iFormat) ) ] );
            end
        end
    end
end

%% Create config file

if createConfigBool
    createConfigFile( currentPath, folder, fileName, commentString, ... % File settings
        plotTimes, plotVeloDistBool, smoothPlotBool, plotDensityBool, plotLogBool, plot2DBool, saveFiguresBool, saveFormat, ... % Preferences
        collisionBool, nMin, kMax, negativeVeloBool, veloSmoothNonColBool, veloSmoothColBool, smoothWidth, keepParticleBool, ... % Computational restrictions
        time, angle, radius, velo, ... % Dimensional limits
        uc, ucRatio, atomUC, nAtomUC, nOxidePerElement, targetDensity, densityRatio, absorbedRatio, heatTarget, energyExcitation, ... % Material parameters
        bg, bgPressure, bgTemperature, bgDensity, ... % Background gas parameters
        spotWidth, spotHeight, ablationDepth, ablationVolume, laserFluence, energyLaser, ... % Laser parameters
        cosPowerFit, distributionType, initVeloDistWidth, ... % Model parameters
        nAtomAblated, nPlasmaTotal, nBGTotal, nNonColTotal, nColTotal, ... % Number of particles
        plasmaTempStart, plasmaTempEnd, bgTempStart, bgTempEnd, ... % Temperatures
        debugBool, durationTotal, nLoops ); % Technical                         
end
