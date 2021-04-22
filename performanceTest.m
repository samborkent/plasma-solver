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

tic;

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
nMin = 1000;

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 8E-6;     % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 10E-5;     % Radial step size [m]

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
bgPressure      = 0.02E2;

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
%   * Calculate per target for most accurate result.
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

%   Metal atom              M   : 1
%   Oxygen atom             O   : 2

nSpecies = nElements + 2*(nElements - 1);

% Plasma matrices
plasmaMatrix = zeros(nVelo, nSpecies, nRadius);
plasmaSub = plasmaMatrix;
plasmaAdd = plasmaSub;

% Background matrices
bgMatrix = zeros(nVelo, nRadius);
bgSub = bgMatrix;
bgAdd = bgSub;

%% Pre-calculations

% Matrix containing all possible velocity weights
veloWeight = (velo' - velo) ./ (velo' + velo);

% Angular plasma particle distribution
nParticleAngle = angularDistribution( angle,        ...
                                      radiusDelta,  ...
                                      cosPowerFit,  ...
                                      nAtomAblated );

%% Initialize figures

% % 1D plasma propagation figure
% figPropagation1D = figure('Name', '1D plasma particle propagation');
% hold on;
% 
% % 1D background propagation figure
% figBGPropagation1D = figure('Name', '1D background particle propagation');
% hold on;

%% Main program

% Temporary values for testing
radiusBG = 2*PT.O.RADIUS;
sigmaArray   = pi .* ( radiusBG + [uc.ELEMENTS.RADIUS] ).^2;
massBG    = 2*PT.O.MASS;
massArray = [uc.ELEMENTS.MASS];

iNewVeloMatrix = zeros(nVelo, nElements, nVelo);
iNewVeloBGMatrix = zeros(nVelo, nElements, nVelo);

for iSpecies = 1 : nElements
    iNewVeloMatrix(:, iSpecies, :) = round( (velo'.*(massArray(iSpecies) - massBG) ...
        + 2*massBG.*velo) ./ (massArray(iSpecies) + massBG) ./ veloDelta );
    
    iNewVeloBGMatrix(:, iSpecies, :) = round( (velo'.*(massBG - massArray(iSpecies)) ...
        + 2*massArray(iSpecies).*velo) ./ (massArray(iSpecies) + massBG) ./ veloDelta );
end

% Prevent index out-of-range error
iNewVeloMatrix(iNewVeloMatrix < 1) = 1;
iNewVeloBGMatrix(iNewVeloBGMatrix < 1) = 1;
iNewVeloMatrix(iNewVeloMatrix > nVelo) = nVelo;
iNewVeloBGMatrix(iNewVeloBGMatrix > nVelo) = nVelo;

nMinSum = nMin * nElements;
                                  
for iAngle = 1 : 1 % nAngle
%% Calculations per angle

% Skip angle if there are not enough particles present
if nParticleAngle(iAngle) < nMin
    continue
end

% Reset matrices
if iAngle > 1
    plasmaMatrix = plasmaMatrix.*0;
    bgMatrix = bgMatrix.*0;
end

% Pre-allocate background gas particle matrix and fill matrix based on
%   calculated density, also retrieve bin volumes per radial bin
[bgMatrix(1,:), binVolume] = fillBGMatrix( bgDensity,    ...
                                      radius,       ...
                                      angle,        ...
                                      iAngle,       ...
                                      radiusDelta );
                                  
% Total number of background particles at this angle
nBGTotal = sum( bgMatrix(1, :) );

% Total number of plasma particles at this angle
nPlasmaTotal = nParticleAngle(iAngle) .* uc.AMOUNT ./ sum(uc.AMOUNT);
                                  
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
                                           
% Fill into the plasma matrix
plasmaMatrix(:, 1:nElements, 1) = nParticleVeloInit;

for iTime = 1 : nTime
%% Calculations per time step

%--------------------------------------------------------------------------
% Reset update matrices for t+1
%--------------------------------------------------------------------------

plasmaSub = plasmaSub.*0;
plasmaAdd = plasmaSub;
bgSub = bgSub.*0;
bgAdd = bgSub;

%--------------------------------------------------------------------------
% Calculate furthest possible occupied radial bin
%--------------------------------------------------------------------------

% Only loop through radial bins where particles could potentially be
%   present
if (nVelo * iTime) < nRadius
    endRadius = (nVelo * iTime) - 1;
else
    endRadius = nRadius - 1;
end

for iSpecies = nElements : -1 : 1
% Loop from lowest to highest mass

for iRadius = endRadius : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

% Calculate sum of all particles in radial bin
nPlasmaRadius = sum( plasmaMatrix(:, iSpecies, iRadius) );

% If no particles are present in radial bin, skip to next radial bin
if nPlasmaRadius < nMin
    continue
end

for iVelo = nVelo : -1 : 2
%% Calculations per plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Only loop through velocities that traverse at least one radial bin per
%   time step

%--------------------------------------------------------------------------
% Get plasma particle variables
%--------------------------------------------------------------------------

% Number of plasma particles in current bin
nPlasma = plasmaMatrix(iVelo, iSpecies, iRadius);

% If no filled bins remain, skip to next radial bin
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
% Calculate number of collisions
%--------------------------------------------------------------------------

% Background particle density in current bin
bgDensity = nBG / binVolume(thisRadius);

% Collision rate
colRate = ( bgDensity * radiusDelta * veloWeight(iVelo, iVeloBG) ) ...
          * sigmaArray(iSpecies);

% Number of collided particles
nCol = colRate .* nPlasmaTemp;

% Skip iteration if the number of collisions is smaller than the threshold
if nCol < nMin
    continue
end

% Set maximum number of collisions
if nPlasmaTemp > nBG
    % Limit number of collisions to the maximum number of collisions
    nCol(nCol > nBG) = nBG;
else
    % Limit number of collisions to the maximum number of collisions
    nCol(nCol > nPlasmaTemp) = nPlasmaTemp;
end

% Update number of non-collided plasma particles
nPlasmaTemp = nPlasmaTemp - nCol;

%--------------------------------------------------------------------------
% Calculate new radial positions after collision
%--------------------------------------------------------------------------

% Reshape indice arrays
% iNewVelo = iNewVeloMatrix(iVelo, iSpecies, iVeloBG);
% iNewVeloBG = iNewVeloBGMatrix(iVelo, iSpecies, iVeloBG);

% deltaVelo = iVelo - iVeloBG;
% massSum = massArray(iSpecies) + massBG;
% 
% iNewVelo = iVelo + round( deltaVelo * (massArray(iSpecies) - massBG) / massSum );
% iNewVeloBG = iVeloBG + round( 2 * massArray(iSpecies) * deltaVelo / massSum );
% 
% iNewVelo(iNewVelo < 1) = 1;
% iNewVelo(iNewVelo > nVelo) = nVelo;
% iNewVeloBG(iNewVeloBG < 1) = 1;
% iNewVeloBG(iNewVeloBG > nVelo) = nVelo;

if iVeloBG == 1
    iNewVelo = iVelo + round( iVelo * (massArray(iSpecies) - massBG) / (massArray(iSpecies) + massBG) );
    iNewVeloBG = iVeloBG + round( 2 * massArray(iSpecies) * iVelo / (massArray(iSpecies) + massBG) );
else
    iNewVelo = iVelo + round( (iVelo - iVeloBG) * (massArray(iSpecies) - massBG) / (massArray(iSpecies) + massBG) );
    iNewVeloBG = iVeloBG + round( 2 * massArray(iSpecies) * (iVelo - iVeloBG) / (massArray(iSpecies) + massBG) );
end

if iNewVelo < 1
    iNewVelo = 1;
end
if iNewVelo > nVelo
    iNewVelo = nVelo;
end
if iNewVeloBG < 1
    iNewVeloBG = 1;
end
if iNewVeloBG > nVelo
    iNewVeloBG = nVelo;
end

% New radial index of collided plasma particles
nNewRadiusDelta = round( velo(iNewVelo) .* (nRadiusDelta - jRadius) ...
                         .* timeDelta / (nRadiusDelta * radiusDelta) );
                       
% New radial index of collided background particles
nNewRadiusDeltaBG = round( velo(iNewVeloBG) .* (nRadiusDelta - jRadius) ...
                           .* timeDelta / (nRadiusDelta * radiusDelta) );

% Prevent index out-of-range error
nNewRadiusDelta(thisRadius + nNewRadiusDelta > nRadius) = nRadius - thisRadius;
nNewRadiusDeltaBG(thisRadius + nNewRadiusDeltaBG > nRadius) = nRadius - thisRadius;

%--------------------------------------------------------------------------
% Update matrices
%--------------------------------------------------------------------------

% Remove collided background particles from collision position
bgSub( iVeloBG, thisRadius ) = bgSub( iVeloBG, thisRadius ) - nCol;

% Add collided background particles to new position and new velocity after
%   collision
bgAdd( iNewVeloBG, thisRadius + nNewRadiusDeltaBG ) = ...
    bgAdd( iNewVeloBG, thisRadius + nNewRadiusDeltaBG ) + nCol;

% Remove collided plasma particles from starting position
plasmaSub( iVelo, iSpecies, iRadius ) = ...
    plasmaSub( iVelo, iSpecies, iRadius ) - nCol;

% Add collided plasma particles to new position and new velocity after
%   collision
plasmaAdd( iNewVelo, iSpecies, thisRadius + nNewRadiusDelta ) = ...
   plasmaAdd( iNewVelo, iSpecies, thisRadius + nNewRadiusDelta ) + nCol;

end % Background velocity loop

end % Traversed radial bin loop

end % Non-oxidizing plasma velocity loop

end % Radius loop

end % Species loop

%--------------------------------------------------------------------------
% Update non-colliding particles
%--------------------------------------------------------------------------

% Remove previously moved collided particles
plasmaMatrix = plasmaMatrix + plasmaSub;
bgMatrix = bgMatrix + bgSub;

% Update non-collided plasma particles
plasmaMatrix = updateMatrix( plasmaMatrix, ...
                             nRadius, ...
                             nVelo, ...
                             nMin );

% Update non-collided background particles
bgMatrix = updateMatrix( bgMatrix, ...
                         nRadius, ...
                         nVelo, ...
                         nMin );

% Add back the previously moved collided particles
plasmaMatrix = plasmaMatrix + plasmaAdd;
bgMatrix = bgMatrix + bgAdd;

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

        nBGRadius = nParticlesPerRadius( radius, bgMatrix );

        figure(5);
        hold on;

        plot( radius, nBGRadius ./ binVolume, ...
              'LineWidth', 2, ...
              'DisplayName', [num2str(time(iTime), 3) ' s' ] );

        xlim([0 0.05]);
        legend;
        hold off;
        
        species = [1 2 5];

        % Loop through plasma species
        for iSpecies = 1 : 2
            % Calculate total number of particles per radial bin
            nParticleRadius = nParticlesPerRadius( radius, plasmaMatrix(:, species(iSpecies), :) );
            
%             nParticleRadius = smoothdata(nParticleRadius, 2, 'gaussian', 50);
%             
%             if species(iSpecies) == nSpecies
%                 nParticleNorm = nBGTotal / sum(nParticleRadius);
%             else
%                 nParticleNorm = nPlasmaTotal(iSpecies) / sum(nParticleRadius);
%             end
%             
%             nParticleRadius = nParticleRadius .* nParticleNorm;
            
            figure(iSpecies);
            hold on;
              
            plot( radius, nParticleRadius, ...
                  'LineWidth', 2, ...
                  'DisplayName', [num2str(time(iTime), 3) ' s' ] );
            
            xlim([0 0.05]);
            legend;
            hold off;
        end
        
%         % Smooth data
%         nPlasmaRadius = smoothdata(nPlasmaRadius, 2, 'gaussian', 100);
%         nBGRadius = smoothdata(nBGRadius, 2, 'gaussian', 50);
%                 
%         % Normalization factor to conserve number of particles
%         nPlasmaNorm = nPlasmaTotal / sum(nPlasmaRadius);
%         nBGNorm = nBGTotal / sum(nBGRadius);
%         
%         % Normalize
%         nPlasmaRadius = nPlasmaRadius .* nPlasmaNorm;
%         nBGRadius = nBGRadius .* nBGNorm;

%         % Plot 1D plasma propagation
%         figure(figPropagation1D);
%         plot( radius, nPlasmaRadius, ...
%               'LineWidth', 2, ...
%               'DisplayName', [num2str(time(iTime), 3) ' s' ] );
%           
%         % Plot 1D background propagation
%         figure(figBGPropagation1D);
%         plot( radius, nBGRadius ./ binVolume, ...
%               'LineWidth', 2, ...
%               'DisplayName', [num2str(time(iTime), 3) ' s' ] );
    end % If Time
end % If Angle

end % Time loop

end % Angle loop

% %% Plot settings
% 
% % 1D plasma propagation
% figure(figPropagation1D);
% legend;
% xlim([0 0.05]);
% % ylim([0 5E12]);
% title('1D propagation of Ti from TiO_2 target in 0.02 mbar O_2 background');
% xlabel('Target distance [m]');
% ylabel('Number of particles per radial bin');
% 
% % 1D background propagation
% figure(figBGPropagation1D);
% legend;
% xlim([0 0.05]);
% % ylim([0 1E21]);
% title('1D propagation of background O_2 gas particles');
% xlabel('Target distance [m]');
% ylabel('Particle density [m^-^3]');

toc;
