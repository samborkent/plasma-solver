%% Initialization
% Do not change!

%--------------------------------------------------------------------------
% Clear workspace
%--------------------------------------------------------------------------
clc, clear, close all

%--------------------------------------------------------------------------
% Timer
%--------------------------------------------------------------------------

% Time the entire script
durationTotalTimer = tic;

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
plotVelDisInit = false;

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
nMin = 1E8;

% Maximum number of collisions
kMax = 2;

% Switch velocity smoothing for non -collided particles on or off (true / false)
%   * Circumvents quantization error, impacts performance greatly
%   * Conserves number of particles
veloSmoothingBool = true;

% Width of the Gaussian smoothing function for smoothing the non-collided
%   particle verlocities
%   * 3 corresponds roughly due to 400 m/s
smoothWidth = 3;

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 10E-6;    % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 6E-5;     % Radial step size [m]

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
bgPressure      = 0.1E2; % 0.02E2;

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

nSpecies = nElements + 2*(nElements - 1) + 1;

% Plasma particle matrix
timer0 = tic;
% plasmaMatrix = zeros(kMax, nSpecies, nVelo, nRadius);
particleMatrix = zeros(nSpecies, kMax, nVelo, nRadius);
% plasmaSub = plasmaMatrix.*0;
collisionMatrix = zeros(nSpecies, 1, nVelo, nRadius);
collisionRate = zeros(nSpecies-1, 1, nVelo, nRadius);
time0 = toc(timer0);

%% Pre-calculations

sigma = zeros(nSpecies-1, 1);
mass = zeros(nSpecies-1, 1);

sigma(1:4)   = pi .* ( 2*PT.O.RADIUS + [PT.Ti.RADIUS; PT.O.RADIUS; PT.Ti.RADIUS+PT.O.RADIUS; PT.Ti.RADIUS+2*PT.O.RADIUS] ).^2;
mass(1:4)    = [PT.Ti.MASS PT.O.MASS PT.Ti.MASS+PT.O.MASS PT.Ti.MASS+2*PT.O.MASS];
massBG  = 2*PT.O.MASS;

% Matrix containing all possible velocity weights
veloWeight = (velo' - velo) ./ (velo' + velo);

iNewVeloMatrix = zeros(nVelo, nElements, nVelo);
iNewVeloBGMatrix = zeros(nVelo, nElements, nVelo);

for iSpecies = 1 : 4 % nSpecies-1
    iNewVeloMatrix(:, iSpecies, :) = round( (velo'.*(mass(iSpecies) - massBG) ...
        + 2*massBG.*velo) ./ (mass(iSpecies) + massBG) ./ veloDelta );
    
    iNewVeloBGMatrix(:, iSpecies, :) = round( (velo'.*(massBG - mass(iSpecies)) ...
        + 2*mass(iSpecies).*velo) ./ (mass(iSpecies) + massBG) ./ veloDelta );
end

% Prevent index out-of-range error
iNewVeloMatrix(iNewVeloMatrix < 1) = 1;
iNewVeloBGMatrix(iNewVeloBGMatrix < 1) = 1;
iNewVeloMatrix(iNewVeloMatrix > nVelo) = nVelo;
iNewVeloBGMatrix(iNewVeloBGMatrix > nVelo) = nVelo;

% Angular plasma particle distribution
nParticleAngle = angularDistribution( angle,        ...
                                      radiusDelta,  ...
                                      cosPowerFit,  ...
                                      nAtomAblated );

%% Main program

% Pre-allocate background gas particle matrix and fill matrix based on
%   calculated density, also retrieve bin volumes per radial bin
[particleMatrix(1, 1, 1, :), binVolume] = fillBGMatrix( bgDensity,    ...
                                                     radius,       ...
                                                     angle,        ...
                                                     1,            ...
                                                     radiusDelta );
                                        
% % Total number of background particles at this angle
% nBGTotal = sum( particleMatrix(1, 1, 1, :) );
% 
% % Total number of plasma particles at this angle
% nPlasmaTotal = nParticleAngle(1) ./ [uc.AMOUNT];

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
                                               nParticleAngle(1),       ...
                                               1 );

% Fill into the plasma matrix
particleMatrix(2:3, 1, :, 1) = reshape(nParticleVeloInit(:, 1:2)', nElements, 1, nVelo, 1);

nLoops = 0;
errorCount = 0;

nMinSum = nMin; % nMin * (nSpecies - 1) * kMax;

colorArray = ['b' 'r' 'k'];
plotArray = {'Uncollided', 'Collided', 'Total'};

plotIndex = 0;

for iTime = 1 : 62 % nTime
%% Calculations per time step

% Reset collision matrix
collisionMatrix = collisionMatrix.*0;

% Total number of particles
particleTotal = sum( particleMatrix, 'all' );

% if iTime > 1
% %     % Smooth velocities of non-collided particles
% %     particleMatrix(:, 1, 2:end, :) = smoothdata( particleMatrix(:, 1, 2:end, :), 3, 'gaussian', smoothWidth );
%     % Smooth velocities
%     particleMatrix(:, :, 2:end, :) = smoothdata( particleMatrix(:, :, 2:end, :), 3, 'gaussian', smoothWidth );
% end

% Remove all particles below threshold
particleMatrix(particleMatrix < nMin) = 0;

% Normalize number of particles and add removed particles back to filled bins
particleMatrix = ( particleMatrix ./ sum(particleMatrix, 'all') ) .* particleTotal;

% Calculate collision rate matrix
collisionRate = ( sum(particleMatrix(1, :, :, :), 2) ./ reshape(binVolume, 1, 1, 1, nRadius) ) ...
                .* radiusDelta .* sigma;

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
particleMatrix = alternativeUpdate(particleMatrix, nVelo, nSpecies, kMax);

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
        
        for k = 1 : 2 % kMax
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
        for i = 1 : nRadius
            veloMeanSquared(i) = sum(reshape(sum(particleMatrix(iSpecies, :, :, i), 2), 1, nVelo) .* (velo.^2)) / sum(particleMatrix(iSpecies, :, :, i), 'all');
            
            if iSpecies == 1
                tempMean(i) = (massBG * veloMeanSquared(i)) / (3 * C.BOLTZMANN);
            else
                tempMean(i) = (mass(iSpecies) * veloMeanSquared(i)) / (3 * C.BOLTZMANN);
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

durationTotal = toc(durationTotalTimer);

function matrix = alternativeUpdate(matrix, nVelo, nSpecies, kMax)
    for iVelo = nVelo : -1 : 2
%         matrix(:, :, iVelo, end-iVelo+1) = ...
%             sum( matrix(:, :, iVelo, end-iVelo+1:end), 'all' );
        matrix(:, :, iVelo, iVelo:end) = matrix(:, :, iVelo, 1:end-iVelo+1);
        matrix(:, :, iVelo, 1:iVelo-1) = zeros(nSpecies, kMax, 1, iVelo-1);
    end
end

function matrix = quickUpdate(matrix, nVelo, nRadius, nMin)
    for iRadius = nRadius - 1 : -1 : 1   
        for iVelo = nVelo : -1 : 2
            if matrix(iVelo, iRadius) > nMin
                if (iRadius + iVelo - 1) > nRadius
                    nRadiusTraveled = nRadius - iRadius;
                else
                    nRadiusTraveled = iVelo - 1;
                end
            
                matrix(iVelo, iRadius + nRadiusTraveled) = ...
                    matrix(iVelo, iRadius + nRadiusTraveled) ...
                    + matrix(iVelo, iRadius);
                matrix(iVelo, iRadius) = 0;
            end
        end
    end
end