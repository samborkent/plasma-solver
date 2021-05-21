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

tic;

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
nMin = 1E5;

% Maximum number of collisions
kMax = 3;

% Switch velocity smoothing for non -collided particles on or off (true / false)
%   * Circumvents quantization error, impacts performance greatly
%   * Conserves number of particles
veloSmoothingBool = true;

% Width of the Gaussian smoothing function for smoothing the non-collided
%   particle velocities
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
timeMax     = 4E-6;     % End time [s]
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

colMatrix = zeros(kMax, nVelo, nRadius);

colSub = colMatrix;
colAdd = colMatrix;

% Plasma particle update matrices
plasmaSub = zeros(nVelo, nRadius);  % Removing particles
plasmaAdd = plasmaSub;              % Adding particles

% Background particle matrix
bgMatrix = zeros(nVelo, nRadius);

% Background particle update matrices
bgSub = bgMatrix;   % Removing particles
bgAdd = bgMatrix;   % Adding particles

%% Pre-calculations

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
% sigma   = pi * ( 2*PT.O.RADIUS + PT.Ti.RADIUS )^2;
sigma   = pi * ( 2*PT.O.RADIUS + PT.O.RADIUS )^2;
% mass    = PT.Ti.MASS;
mass    = PT.O.MASS;
massBG  = 2*PT.O.MASS;

index = 0;
                                  
for iAngle = 1 : 1 % nAngle
%% Calculations per angle

% Skip angle if there are not enough particles present
if nParticleAngle(iAngle) < nMin
    continue
end

% Pre-allocate background gas particle matrix and fill matrix based on
%   calculated density, also retrieve bin volumes per radial bin
[bgMatrix(1, :), binVolume] = fillBGMatrix( bgDensity,    ...
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
    colMatrix = colMatrix.*0;
end
                                           
% Fill into the plasma matrix
plasmaMatrix(:, 1) = nParticleVeloInit(:, 2);
colMatrix(1, :, 1) = nParticleVeloInit(:, 2);

for iTime = 1 : nTime
%% Calculations per time step

%--------------------------------------------------------------------------
% Reset update matrices for t+1
%--------------------------------------------------------------------------

plasmaSub = plasmaSub.*0;
plasmaAdd = plasmaSub;
bgSub = bgMatrix.*0;
bgAdd = bgSub;

colSub = colSub.*0;
colAdd = colSub;

for iRadius = (nRadius - 1) : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

for iVelo = nVelo : -1 : 2
%% Calculations per non-oxidizing plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Only loop through velocities that traverse at least one radial bin per
%   time step

%--------------------------------------------------------------------------
% Get plasma particle variables
%--------------------------------------------------------------------------

colWeights = colMatrix(:, iVelo, iRadius);

% Number of plasma particles in current bin
% nPlasma = plasmaMatrix(iVelo, iRadius);
nPlasma = sum( colWeights );

% Skip bins with number of particles below threshold
if nPlasma < nMin
    continue
end

colWeights = colWeights ./ nPlasma;

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
% Calculate velocity weight factors
%--------------------------------------------------------------------------

% Relative velocity weight factors
veloWeight = (velo(iVelo) - velo(iVeloBG)) ./ (velo(iVelo) + velo(iVeloBG));

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
                       
% Prevent index out-of-range error
if thisRadius + nNewRadiusDelta > nRadius
    nNewRadiusDelta = nRadius - thisRadius;
end
if thisRadius + nNewRadiusDeltaBG > nRadius
    nNewRadiusDeltaBG = nRadius - thisRadius;
end

%--------------------------------------------------------------------------
% Update matrices
%--------------------------------------------------------------------------


% Remove collided particles from starting position
plasmaSub(iVelo, iRadius) = plasmaSub(iVelo, iRadius) - nCol;
bgSub(iVeloBG, thisRadius) = bgSub(iVeloBG, thisRadius) - nCol;

nColArray = colWeights .* nCol;

colSub(:, iVelo, iRadius) = colSub(:, iVelo, iRadius) - nColArray;

nColArray = [0;nColArray];
nColArray(end-1) = nColArray(end-1) + nColArray(end);
nColArray(end) = [];

% Add collided particles to new position and new velocity after collision
plasmaAdd(iNewVelo, thisRadius + nNewRadiusDelta) = ...
    plasmaAdd(iNewVelo, thisRadius + nNewRadiusDelta) + nCol;
bgAdd(iNewVeloBG, thisRadius + nNewRadiusDeltaBG) = ...
    bgAdd(iNewVeloBG, thisRadius + nNewRadiusDeltaBG) + nCol;

colAdd(:, iNewVelo, thisRadius + nNewRadiusDelta) = ...
    colAdd(:, iNewVelo, thisRadius + nNewRadiusDelta) + nColArray;

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

colMatrix = colMatrix + colSub;

% Update non-collided particles
plasmaMatrix = updateMatrix( plasmaMatrix, nRadius, nVelo, nMin );
bgMatrix = updateMatrix( bgMatrix, nRadius, nVelo, nMin );

colMatrix = updateCol( colMatrix, nRadius, nVelo, nMin );

% Add back the previously moved collided particles
plasmaMatrix = plasmaMatrix + plasmaAdd;
bgMatrix = bgMatrix + bgAdd;

colMatrix = colMatrix + colAdd;

%--------------------------------------------------------------------------
% Velocity smoothing
%--------------------------------------------------------------------------

% If velocity smoothing is enabled
if veloSmoothingBool
    % Loop through number of collisions
    for k = 1 : kMax
        % Get matrix per number of collisions
        kMatrix = squeeze( colMatrix(k, :, :) );

        % Sum all elements in number of collision matrix
        kMatrixSum = sum( kMatrix, 'all' );

        % If the number of collision matrix is filled
        if kMatrixSum >= nMin
            % Set all bins with number of particles below threshold to zero
            kMatrix(kMatrix < nMin) = 0;

            % For the non-collided particles
            if k == 1
                kMatrix = smoothdata( kMatrix, 1, 'gaussian', smoothWidth );
            end

            % Normalize and multiply by original number of particles before
            %   emptying bins below threshold
            kMatrix = kMatrix .* ( kMatrixSum / sum(kMatrix, 'all') );

            % Insert into full matrix
            colMatrix(k, :, :) = kMatrix;
        end
    end
end

%--------------------------------------------------------------------------
% Plot 1D propagation
%--------------------------------------------------------------------------

% Only for first angle (center of the plume)
if iAngle == 1   
    % Only for specific times
    if (iTime ==  6) || (iTime == 11) || (iTime == 16) || (iTime == 31) || (iTime ==  46) || ...
       (iTime == 51) || (iTime == 61) || (iTime == 71)
%     if (iTime == 11)  || (iTime ==  21) || (iTime == 31) || (iTime ==  41) || ...
%        (iTime == 51)  || (iTime ==  61) || (iTime == 91) || (iTime == 121) || ...
%        (iTime == 151) || (iTime == 181)

        % Increment subplot index
        index = index + 1;

        % Initialize sum of all collisions array
        nPlasmaRadiusSum = zeros(1, nRadius);
        
        for k = 1 : kMax
            % Calculate number of particles with this number of collisions
            %   per radial bin
            nPlasmaRadius = nParticlesPerRadius( radius, colMatrix(k, :, :) );
            
            if ~veloSmoothingBool
                if k == 1
                    % Total number of particles for this number of collisions
                    nPlasmaPerK = sum( colMatrix(k, :, :), 'all' );

                    % Fit a normalized Gaussian curve to the number of particles
                    nPlasmaRadius = fitGaussian( radius, nPlasmaRadius, nMin );

                    % Multiply by number of particles
                    nPlasmaRadius = nPlasmaRadius .* nPlasmaPerK;
                else
                    % Smooth data
                    nPlasmaRadius = smoothdata( nPlasmaRadius, 2, 'gaussian', 50 );
                end
            end
            
            % Add to sum of all collisions array
            nPlasmaRadiusSum = nPlasmaRadiusSum + nPlasmaRadius;
            
            % Plot 1D plasma propagation
            figure(figPropagation1D);
            subplot(1, 4, index);
            hold on;
            plot( radius, nPlasmaRadius, ...
                  'LineWidth', 2, ...
                  'DisplayName', [num2str(k - 1) ' coll.' ] );
        end
        
        plot( radius, nPlasmaRadiusSum, ...
          'LineWidth', 2, ...
          'Color', 'k', ...
          'DisplayName', 'Sum' );
      
        title([num2str(time(iTime), 3) ' s' ]);
        xlim([0 0.05]);
        ylim([0 24E11]);
        xlabel('Distance [m]');
        if index == 1
            ylabel('Number of particles');
        end

        % Calculate total number of background particles per radial bin
        nBGRadius = nParticlesPerRadius( radius, bgMatrix );
        
        % Smooth background data
        nBGRadius = smoothdata( nBGRadius, 2, 'gaussian', 50 );
          
        % Plot 1D background propagation
        figure(figBGPropagation1D);
        plot( radius, nBGRadius ./ binVolume, ...
              'LineWidth', 2, ...
              'DisplayName', [num2str(time(iTime), 3) ' s' ] );
        
    end
end

%--------------------------------------------------------------------------
% Debug per time step
%--------------------------------------------------------------------------

if debugBool
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
% sgtitle('1D propagation of Ti from TiO_2 target in 0.02 mbar O_2 background');
sgtitle('1D propagation of O from TiO_2 target in 0.02 mbar O_2 background');

% 1D background propagation
figure(figBGPropagation1D);
legend;
xlim([0 0.05]);
% ylim([0 1E21]);
% title('1D propagation of background O_2 gas particles from collisions with Ti plasma at 0.02 mbar');
title('1D propagation of background O_2 gas particles from collisions with O plasma at 0.02 mbar');
xlabel('Target distance [m]');
ylabel('Particle density [m^-^3]');

toc;
