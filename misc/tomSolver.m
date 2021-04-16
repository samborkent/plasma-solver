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

% Maximum number of average collisions per bin
colMax = 3;

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

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 12E-6;     % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 5E-5;     % Radial step size [m]

% Velocity limits
veloMax         = 2.5E4; % Maximal initial velocity
veloStepsize    = 1;     % Velocity step size

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
angle   = angleMin                  : angleDelta    : angleMax;
nAngle  = numel(angle);

% Temporal axis
time    = timeMin                   : timeDelta     : timeMax - timeDelta;
nTime   = numel(time);

% Radial axis
radius  = radiusMin + radiusDelta   : radiusDelta   : radiusMax;
nRadius = numel(radius);

% Velocity array
veloDelta   = (radiusDelta / timeDelta) / veloStepsize;
velo        = 0 : veloDelta : veloMax - veloDelta;
nVelo       = numel(velo);

% Number of radial bins traveled in one time step per velocity bin
nRadiusTraveled = round( (velo .* timeDelta) ./ radiusDelta );

% Index of first velocity high enough to traverse one radial bin per
%   timestep
iFirstVelo = find(nRadiusTraveled, 1);

%% Pre-allocations

% Number field of the plasma matric
numFields = 2;

% Indexing of the plasma matrix fields
nField = 1;     % Number of particles
kField = 2;     % Number of collisions

% Plasma particle matrix
plasmaMatrix = zeros(nVelo, nRadius, numFields);

% Collision array
tempRad = zeros(1, nRadius);

bgNew = zeros(nVelo, nRadius);
plasmaNew = zeros(nVelo, nRadius);

%% Initialize figures

% 1D plasma propagation figure
figPropagation1D = figure('Name', '1D Plasma propagation');
hold on;

% 1D background propagation figure
figBGPropagation1D = figure('Name', '1D Background propagation');
hold on;

% 1D plasma propagation with separated collisions
figCollisions1D = figure('Name', '1D Plasma propagation: seperated number of collisions');

% Index for number of plotted times
plotIndex = 0;

%% Angular plasma particle distribution

nParticleAngle = angularDistribution( angle,        ...
                                      radiusDelta,  ...
                                      cosPowerFit,  ...
                                      nAtomAblated );

%% Main program

% Temporary values for testing
sigma   = pi * ( 2*PT.O.RADIUS + PT.Ti.RADIUS )^2;
mass    = PT.Ti.MASS;
massBG  = 2*PT.O.MASS;
                                  
for iAngle = 1 : 1 % (nAngle - 1)
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
                                      iAngle );

% Calculate the initial particle velocity distribution
[nPlasmaVeloInit, veloPlasma] = initialVelocityDistribution( false,     ...
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
                                           
nPlasmaTotal = sum(nPlasmaVeloInit(:, 1));
                                           
iVeloPlasma = round(veloPlasma * timeDelta / radiusDelta);

% Fill into the plasma matrix
% plasmaMatrix(:, 1, nField) = nPlasmaVeloInit(:, 1);
plasmaMatrix(iVeloPlasma, 1, nField) = sum(nPlasmaVeloInit(:, 1));

% [DEBUG] Check initial total amount particles in matrix 
startPlasma = sum(sum(plasmaMatrix(:,:,nField)));
startBG = sum(sum(bgMatrix));

for iTime = 1 : nTime
%% Calculations per time step

% Plasma matrix at t+1
plasmaNew = zeros(nVelo, nRadius, numFields);

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
nPlasma = plasmaMatrix(iVelo, iRadius, nField);

if nPlasma < nMin
    continue
end

nPlasmaTemp = nPlasma;

%--------------------------------------------------------------------------
% Set number of traveled radial bins
%--------------------------------------------------------------------------

% Restrict number of traveled bins to the total number of
%   radial bins
if ( iRadius + nRadiusTraveled(iVelo) ) > nRadius
    nRadiusDelta = nRadius - iRadius;
else
    nRadiusDelta = nRadiusTraveled(iVelo);
end

for thisRadius = 0 : (nRadiusDelta - 1)
%% Calculations per traversed radial bin
% Loop from current bin to one bin before final bin

    for jVelo = 1 : (iVelo - 1)
        %% Calculations per background velocity bin
        % Only include filled background velocities smaller than the
        %   plasma velocity
        
        nBGVelo = bgMatrix(jVelo, iRadius + thisRadius);
        
        if nBGVelo < nMin
            continue
        end

        % Background particle density in current bin
        bgDensity = nBGVelo(jVelo) / binVolume(iRadius + thisRadius);

        % Collision rate
        colRate = bgDensity * sigma * radiusDelta * nPlasmaTemp;

        % Number of collided particles
        nCol = colRate * nPlasmaTemp;
        
        if nPlasmaTemp > nBGVelo(jVelo)
            nColMax = nBGVelo(jVelo);
        else
            nColMax = nPlasmaTemp;
        end

        % Limit the number of collided particles to the number of particles
        %   in the bin
        if nCol > nColMax
            nCol = nColMax;
        end
        
        

        %------------------------------------------------------------------
        % Calculate new velocities after collision
        %------------------------------------------------------------------

        if jVelo == 1
            % Calculate new plasma velocity index after plasma-bg collision (eq. 6)
            %   * Only valid if veloPlasma >> veloBG
            iNewVelo = round( iVelo * (mass - massBG) ...
                                    / (mass + massBG) );

            % Calculate new bg velocity index after plasma-bg collision (eq. 7)
            %   * Only valid if veloPlasma >> veloBG
            iNewBGVelo = round( 2 * mass * iVelo ...
                                 / (mass + massBG) );
        else
            iNewVelo = round( (iVelo*(mass - massBG) ...
                               + 2 * massBG * jVelo) ...
                               / (mass + massBG) );
                           
            iNewBGVelo = round( (jVelo*(massBG - mass) ...
                                 + 2 * mass * iVelo) ...
                                 / (mass + massBG) );
        end

        % Prevent index out-of-range error
        if iNewVelo < 1
            iNewVelo = 1;
        elseif iNewVelo > nVelo
            iNewVelo = nVelo;
        end
        if iNewBGVelo < 1
            iNewBGVelo = 1;
        elseif iNewBGVelo > nVelo
            iNewBGVelo = nVelo;
        end

        %----------------------------------------------------------------------
        % Calculate new radial positions after collision
        %----------------------------------------------------------------------

        % New radial indices of collided particles
        iNewPlasmaRadius = round( (iVelo * iThisRadius ...
                           + iNewVelo * (nRadiusDelta - iThisRadius) ) ...
                           * timeDelta / (nRadiusDelta * radiusDelta) );

        if nBGRadiusDelta > 0
            iNewBGRadius(iBGVelo, iThisRadius) = ...
                round( (iBGVelo * iThisRadius ...
                        + iNewVelo * (nBGRadiusDelta - iThisRadius) ) ...
                        * timeDelta / (nBGRadiusDelta * radiusDelta) );
        else
            iNewBGRadius(iBGVelo, iThisRadius) = ...
                round( (iBGVelo * iThisRadius ...
                        + iNewVelo * (nBGRadiusDelta - iThisRadius) ) ...
                        * timeDelta / (nBGRadiusDelta * radiusDelta) );
        end

        % Prevent index out-of-range error
        if iRadius + iNewPlasmaRadius(iBGVelo, iThisRadius) > nRadius
            iNewPlasmaRadius(iBGVelo, iThisRadius) = nRadius - iRadius;
        end
        if iRadius + iNewBGRadius(iBGVelo, iThisRadius) > nRadius
            iNewBGRadius(iBGVelo, iThisRadius) = nRadius - iRadius;
        end

    end % Background velocity loop

    % Total number of collided particles in this radial bin
    nColSum(iThisRadius) = sum(nCol(:, iThisRadius));

    % If the total number of collided particles is greater than the number
    %   of plasma particles
    if nColSum(iThisRadius) > nPlasmaTemp(iThisRadius)

        % Loop through the collision array
        for jVelo = 1 : bgVeloNumel
            % Background velocity index
            iBGVelo = bgVeloArray(jVelo);
            
            % Normalize the number of collisions by the number of plasma
            %   particles
            nCol(iBGVelo, iThisRadius) = ( nCol(iBGVelo, iThisRadius) ...
                                         / nColSum(iThisRadius) ) ...
                                         * nPlasmaTemp(iThisRadius);
        end

        % Set the total number of collisions equal to the number of plasma
        %   particles
        nColSum(iThisRadius) = nPlasmaTemp(iThisRadius);
    end

    % Move non-collided particles to next traveled bin
    nPlasmaTemp(iThisRadius + 1) = nPlasmaTemp(iThisRadius) - nColSum(iThisRadius);

end % Traveled radial bin loop

% Total number of non-colliding particles
nNonColPlasma = nPlasma - sum(nColSum);

% New radial index of uncollided plasma particles
iNewPlasmaRadiusNonCol = iRadius + nRadiusDelta;

% Prevent index out-of-range error
if iNewPlasmaRadiusNonCol > nRadius
    iNewPlasmaRadiusNonCol = nRadius;
end

%--------------------------------------------------------------------------
% Update non-colliding particles
%--------------------------------------------------------------------------

% Add non-colliding plasma particles to new position in t+1 matrix
plasmaNew(iVelo, iNewPlasmaRadiusNonCol, nField) = nNonColPlasma;

for thisRadius = iRadius : iRadius + nRadiusDelta - 1
%% Calculations per traversed radial bin
% Loop from current bin to one bin before final bin

    iThisRadius = thisRadius - iRadius + 1;

    %--------------------------------------------------------------------------
    % Get filled background particle bins
    %--------------------------------------------------------------------------

    % Find velocity bins containing particles below the plasma velocity
    bgVeloArray = find( bgMatrix(1:iVelo-1, thisRadius) );

    % Skip to next radial bin if none of the velocity bins are filled
    if isempty(bgVeloArray)
        continue
    end

    % Corresponding number of background particles
    nBGVelo = bgMatrix(bgVeloArray, thisRadius);

    % Only velocity bins with number of particles above threshold
    bgVeloArray = bgVeloArray(nBGVelo >= nMin);

    % Skip to next radial bin if none of the velocity bins has enough particles
    if isempty(bgVeloArray)
        continue
    end

    % Update number of background particles array
    nBGVelo = nBGVelo(nBGVelo >= nMin);

    for jVelo = 1 : numel(bgVeloArray)
        %% Calculations per background velocity bin
        % Only include filled background velocities smaller than the
        %   plasma velocity

        % Background velocity index
        iBGVelo = bgVeloArray(jVelo);

        %----------------------------------------------------------------------
        % Update colliding plasma particles
        %----------------------------------------------------------------------

        plasmaMatrix(iVelo, iRadius, nField) = ...
            plasmaMatrix(iVelo, iRadius, nField) ...
            - nCol(iBGVelo, iThisRadius);
        
        plasmaMatrix(iNewVelo, iRadius + iNewPlasmaRadius(iBGVelo, iThisRadius), nField) = ...
            plasmaMatrix(iNewVelo, iRadius + iNewPlasmaRadius(iBGVelo, iThisRadius), nField) ...
            + nCol(iBGVelo, iThisRadius);
        
        plasmaMatrix(iNewVelo, iRadius + iNewPlasmaRadius(jVelo, iThisRadius), kField) = 1 ...
            + plasmaNew(iNewVelo, iRadius + iNewPlasmaRadius(jVelo, iThisRadius), kField);
        
%         plasmaNew(iNewVelo, iRadius + iNewPlasmaRadius(jVelo, iThisRadius), nField) = ...
%             nCol(jVelo, iThisRadius);
% 
%         plasmaNew(iNewVelo, iRadius + iNewPlasmaRadius(jVelo, iThisRadius), kField) = 1 ...
%             + plasmaNew(iNewVelo, iRadius + iNewPlasmaRadius(jVelo, iThisRadius), kField);

        %----------------------------------------------------------------------
        % Update colliding background particles velocities
        %----------------------------------------------------------------------
        
        bgMatrix(iBGVelo, thisRadius) = bgMatrix(iBGVelo, thisRadius) ...
            - nCol(iBGVelo, iThisRadius);

        bgMatrix(iNewBGVelo, iRadius + iNewBGRadius(iBGVelo, iThisRadius)) = ...
            bgMatrix(iNewBGVelo, iRadius + iNewBGRadius(iBGVelo, iThisRadius)) ...
            + nCol(iBGVelo, iThisRadius);

    end % Background velocity loop

end % Traveled radial bin loop

% %--------------------------------------------------------------------------
% % Update non-colliding background particles positions
% %--------------------------------------------------------------------------
% 
% nBG = bgMatrix(iVelo, iRadius);
% 
% if nBG < nMin
%     continue
% end
% 
% % Move background particles to new position
% bgMatrix(iVelo, iRadius + nRadiusDelta) = nBG ...
%     + bgMatrix(iVelo, iRadius + nRadiusDelta);
% 
% % Empty previous position
% bgMatrix(iVelo, iRadius) = bgMatrix(iVelo, iRadius) - nBG;

end % Plasma velocity loop

end % Plasma radial loop

% Increment plasma matrix with one time step
plasmaMatrix = plasmaNew;

 %% Plot 1D propagation
    % Only for first angle (center of the plume)
    if iAngle == 1        
        % Only for specific times
        if (iTime == 0)  || (iTime == 1)  || ...
           (iTime == 5)  || (iTime == 11) || ...
           (iTime == 16) || (iTime == 31) || ...
           (iTime == 61) || (iTime == nTime)

            % Calculate total number of particles per radial bin
            nPlasmaRadius = nPlasmaParticlesPerRadius( radius, plasmaMatrix );
            nBGRadius = nPlasmaParticlesPerRadius( radius, bgMatrix );

%             plasmaFit = fit(radius, nPlasmaRadius, 'poly1');
            
            % Normalization factor
            plasmaNorm = startPlasma / sum(nPlasmaRadius);
            bgNorm = startBG / sum(nBGRadius);

            % Normalize to conserve number of particles
            nPlasmaRadius = nPlasmaRadius .* plasmaNorm;
            nBGRadius = nBGRadius .* bgNorm;

            % Plot 1D plasma propagation
            figure(figPropagation1D);
            plot( radius, nPlasmaRadius, ...
                 'DisplayName', [num2str(time(iTime), 3) ' s'] );
             
            % Plot 1D background propagation
            figure(figBGPropagation1D);
            plot( radius, nBGRadius ./ binVolume, ...
                 'DisplayName', [num2str(time(iTime), 3) ' s'] );
            
%             % Increment plot index
%             plotIndex = plotIndex + 1;
%             
%             % Plot 1D plasma propagation with separated collisions
%             figure(figCollisions1D);
%             subplot(1, 5, plotIndex);
%             ylim([0 1E13]);
%             title([num2str(time(iTime) * 1E6) ' \mus']);
%             hold on;
%             
%             % Collision matrix
%             colMatrix = round( plasmaMatrix(:, :, kField) );
%             
%             % Temporary matrix holding number of particles
%             plasmaTemp = plasmaMatrix(:, :, nField);
%             
%             % Loop through number of collisions
%             for iNCol = 0 : colMax
%                 
%                 % Loop through radial bins
%                 for iRadius = 1 : nRadius
%                     % Get all bins with the current number of collisions
%                     tempMatrix = plasmaTemp(colMatrix(:, iRadius) == iNCol, iRadius);
% 
%                     % Sum all particles with current number of collisions
%                     tempRad(iRadius) = sum(tempMatrix);
%                 end
% 
%                 % Plot the number of particles at current time step with
%                 %   current number of collisions
%                 plot(radius, tempRad);
%             end
%             
%             % Plot legend only for the last plot
%             if plotIndex == 1
%                 ylabel('Number of particles per radial bin');
%             elseif plotIndex == 3
%                 xlabel('Target distance [m]');
%             elseif plotIndex == 5
%                 legend('Uncollided', '1 collision', '2 collisions');
%             end
            
        end
    end 

end % Temporal loop

end % Angular loop

%% Plot settings

% 1D plasma propagation
figure(figPropagation1D);
legend;
% xlim([0 0.05]);
ylim([0 1E13]);
title('1D propagation of Ti from TiO_2 target');
xlabel('Target distance [m]');
ylabel('Number of particles per radial bin');
% ylabel('Particle density [m^-^3]');

% 1D background propagation
figure(figBGPropagation1D);
legend;
xlim([0 0.05]);
% ylim([0 1E13]);
title('1D propagation of background O_2 gas particles');
xlabel('Target distance [m]');
% ylabel('Number of particles per radial bin');
ylabel('Particle density [m^-^3]');

