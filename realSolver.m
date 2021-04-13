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
timeMax     = 16E-6;     % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 2E-5;     % Radial step size [m]

% Velocity limits
veloMax         = 2E4; % Maximal initial velocity
veloStepsize    = 1;     % Velocity step size

%--------------------------------------------------------------------------
% Model parameters
%--------------------------------------------------------------------------

% Fit parameter of the cosines power in the angular plasma particle
%   distibution function
cosPowerFit = 16;

% Initial velocity distribution width
initVeloDisWidth = 1500;

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Unit cell
%   * See UnitCells
%   * Options: TiO2, SrTiO3, Li2O, LiO2, LiTi2O4, Li2TiO3, Li4Ti5O12
%   * Preface with UC.
uc = UC.TiO2;

% Absorption coefficient
%   Default: 0.6
absorption = 0.6;

%--------------------------------------------------------------------------
% Background gas parameters
%--------------------------------------------------------------------------

% Background gas pressure during depostion [Pa]
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

% Plasma particle matrix
plasmaMatrix = zeros(nVelo, nRadius);

% Plasma particle update matrices
plasmaSub = plasmaMatrix;   % Removing particles
plasmaAdd = plasmaMatrix;   % Adding particles

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

%% Main program

% Temporary values for testing
sigma   = pi * ( 2*PT.O.RADIUS + PT.Ti.RADIUS )^2;
mass    = PT.Ti.MASS;
massBG  = 2*PT.O.MASS;
                                  
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
nBGTotal = sum(sum( bgMatrix ));

% Total number of plasma particles at this angle
nPlasmaTotal = nParticleAngle(iAngle) / sum(uc.AMOUNT);
                                  
% Calculate the initial particle velocity distribution
[nParticleVeloInit, veloPlasma] = initialVelocityDistribution( true,    ...
                                               false,                   ...
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

% Index of average velocity
iVeloPlasma = round( veloPlasma / veloDelta ) + 1;

% Fill into the plasma matrix
plasmaMatrix(:, 1) = nParticleVeloInit(:, 1);
% plasmaMatrix(iVeloPlasma, 1) = nPlasmaTotal;

for iTime = 1 : nTime
%% Calculations per time step

% [DEBUG]
nBGTotalNew = sum(sum( bgMatrix ));
nPlasmaTotalNew = sum(sum( plasmaMatrix ));

if (nBGTotal - nBGTotalNew) > nMin
    disp(['Background particles not conserved ' ...
          num2str(nBGTotal - nBGTotalNew, '%.3E') ...
          ' particles lost.']);
elseif (nBGTotal - nBGTotalNew) < -nMin
    disp(['Background particles not conserved ' ...
          num2str(nBGTotalNew - nBGTotal, '%.3E') ...
          ' particles gained.']);
end

if (nPlasmaTotal - nPlasmaTotalNew) > nMin
    disp(['Plasma particles not conserved ' ...
          num2str(nParticleAngle(iAngle) - nPlasmaTotalNew, '%.3E') ...
          ' particles lost.']);
elseif (nPlasmaTotal - nPlasmaTotalNew) < -nMin
    disp(['Plasma particles not conserved ' ...
          num2str(nPlasmaTotalNew - nParticleAngle(iAngle), '%.3E') ...
          ' particles gained.']);
end

%--------------------------------------------------------------------------
% Reset update matrices for t+1
%--------------------------------------------------------------------------

plasmaSub = plasmaSub.*0;
plasmaAdd = plasmaSub;
bgSub = bgMatrix.*0;
bgAdd = bgSub;

for iRadius = (nRadius - 1) : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

for iVelo = nVelo : -1 : iFirstVelo
%% Calculations per plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Only loop through velocities that traverse at least one radial bin per
%   time step

%--------------------------------------------------------------------------
% Get plasma particle variables
%--------------------------------------------------------------------------

% Number of plasma particles in current bin
nPlasma = plasmaMatrix(iVelo, iRadius);

% [DEBUG] Warning for negative number of particles
if nPlasma < 0
    disp(['Negative number of particles at: ( ' ...
          num2str(iTime) ', ' num2str(iRadius) ', ' num2str(iVelo) ' )' ]);
end

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

% Calculate total number of particles in this radial bin
thisNPlasma = sum( plasmaMatrix(:, thisRadius) );   % Plasma
thisNBG = sum( bgMatrix(:, thisRadius) );           % Background

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
nBG = bgMatrix(iVeloBG, thisRadius);

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

if bgDensity < 0
    bgDensity = 0;
end

% Normalization factor for the collision rate:
%   Scale by the fraction of non-colliding plasma particles compared to the
%   total number of plasma particles in this radial bin
colNorm = nPlasmaTemp / (nPlasmaTemp + thisNPlasma);

% Collision rate
colRate = bgDensity * sigma * radiusDelta * colNorm;

% Normalization factor for the number of collisions:
%   Scale by the proportion of background particles of current verlocity
%   compared to all background particles in current radial bin
nColNorm = nBG / thisNBG;

% Number of collided particles
nCol = colRate * nPlasmaTemp * nColNorm;

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

% [DEBUG]
if iNewVelo > iVelo
    disp('Energy gained by heavy mass particle from collision.');
end

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
nNewRadiusDelta = round( (iVelo * jRadius ...
                         + iNewVelo * (nRadiusDelta - jRadius) ) ...
                         * timeDelta / (nRadiusDelta * radiusDelta) );
                       
% New radial index of collided background particles
nNewRadiusDeltaBG = round( (iVeloBG * jRadius ...
                           + iNewVeloBG * (nRadiusDelta - jRadius) ) ...
                           * timeDelta / (nRadiusDelta * radiusDelta) );
                       
% Prevent index out-of-range error
if iRadius + nNewRadiusDelta > nRadius
    nNewRadiusDelta = nRadius - iRadius;
end
if iRadius + nNewRadiusDeltaBG > nRadius
    nNewRadiusDeltaBG = nRadius - iRadius;
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
plasmaAdd(iNewVelo, iRadius + nNewRadiusDelta) = ...
    plasmaAdd(iNewVelo, iRadius + nNewRadiusDelta) + nCol;
bgAdd(iNewVeloBG, iRadius + nNewRadiusDeltaBG) = ...
    bgAdd(iNewVeloBG, iRadius + nNewRadiusDeltaBG) + nCol;

end % Background velocity loop

end % Traversed radial bin loop

end % Plasma velocity loop

end % Radius loop

%--------------------------------------------------------------------------
% Update non-colliding particles
%--------------------------------------------------------------------------

% Remove previously moved collided particles
plasmaMatrix = plasmaMatrix + plasmaSub;
bgMatrix = bgMatrix + bgSub;

% % Update non-collided particles
plasmaMatrix = updateMatrix( plasmaMatrix, nRadius, nVelo, ...
                             iFirstVelo, nMin, nRadiusTraveled );
bgMatrix = updateMatrix( bgMatrix, nRadius, nVelo, ...
                         iFirstVelo, nMin, nRadiusTraveled );

% Add back the previously moved collided particles
plasmaMatrix = plasmaMatrix + plasmaAdd;
bgMatrix = bgMatrix + bgAdd;

%--------------------------------------------------------------------------
% Plot 1D propagation
%--------------------------------------------------------------------------

% Only for first angle (center of the plume)
if iAngle == 1        
    % Only for specific times
    if (iTime == 6) || (iTime == 11) || (iTime == 16) || (iTime == 31)

        % Calculate total number of particles per radial bin
        nPlasmaRadius = nParticlesPerRadius( radius, plasmaMatrix );
        nBGRadius = nParticlesPerRadius( radius, bgMatrix );
        
        % Normalization factor to conserve number of particles
        plasmaNorm = nPlasmaTotal / sum(nPlasmaRadius);
        bgNorm = nBGTotal / sum(nBGRadius);
        
        % Normalize
        nPlasmaRadius = nPlasmaRadius .* plasmaNorm;
        nBGRadius = nBGRadius .* bgNorm;
        
%         % Fit a normalized Gaussian curve to the number of particles
%         nPlasmaFit = fitGaussian( radius, nPlasmaRadius, nMin );
%         
%         % Multiply by the total number of particles at this angle
%         nPlasmaRadius = nPlasmaFit .* nParticleAngle(iAngle);

        % Plot 1D plasma propagation
        figure(figPropagation1D);
        plot( radius, nPlasmaRadius, ...
              'DisplayName', [num2str(time(iTime), 3) ' s'] );
          
        % Plot 1D plasma propagation
        figure(figBGPropagation1D);
        plot( radius, nBGRadius ./ binVolume, ...
              'DisplayName', [num2str(time(iTime), 3) ' s'] );
        
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
title('1D propagation of Ti from TiO_2 target');
xlabel('Target distance [m]');
ylabel('Number of particles per radial bin');

% 1D background propagation
figure(figBGPropagation1D);
legend;
% xlim([0 0.05]);
% ylim([0 1E21]);
title('1D propagation of background O_2 gas particles');
xlabel('Target distance [m]');
ylabel('Particle density [m^-^3]');
