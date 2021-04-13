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
%   bin     : Computational bin with similar position and velocity
%   n       : Number of particles/bins
%   velo    : Velocity
%   uc      : Unit cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Debug settings
%--------------------------------------------------------------------------

% Enable / disable debug modus (true / false)
debugBool = true;

%--------------------------------------------------------------------------
% Plot settings
%--------------------------------------------------------------------------

% Plot initial velocity distribution (true / false)
plotVelDisInit = true;

% Save the initial velocity distribution plot (true / false)
saveVelDisInit = false;

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

numFields = 2;

nField = 1;
kField = 2;

% Plasma particle matrix
plasmaMatrix = zeros(nVelo, nRadius, numFields);

% Collision array
tempRad = zeros(1, nRadius);

%% Angular plasma particle distribution

nParticleAngle = angularDistribution( angle,        ...
                                      radiusDelta,  ...
                                      cosPowerFit,  ...
                                      nAtomAblated );

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

%% DEBUG

if debugBool
    % [DEBUG] Loop counter
    loopCount = 0;

    % [DEBUG]
    nTimesMultipleBGVelo = 0;
end

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
nPlasmaVeloInit = initialVelocityDistribution( plotVelDisInit,          ...
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

% Only plot and save initial particle velocity distribution for the
%   center of the plume
if iAngle == 1
    plotVelDisInit = false;
    saveVelDitInit = false;
end

% Fill into the plasma matrix
plasmaMatrix(:, 1, nField) = nPlasmaVeloInit(:, 1);

% [DEBUG] Check initial total amount particles in matrix 
startPlasma = sum(sum(plasmaMatrix(:,:,nField)));
startBG = sum(sum(bgMatrix));

for iTime = 1 : nTime
    %% Calculations per time step
    
    %----------------------------------------------------------------------
    % Second pass to update background positions
    %----------------------------------------------------------------------
    % Should be unnecessary -> Include in main loop
    
    for iRadius = (nRadius - 1) : -1 : 1
    for iVelo = nVelo : -1 : iFirstVelo
        %------------------------------------------------------------------
        % Set number of traveled radial bins
        %------------------------------------------------------------------

        % Restrict number of traveled bins to the total number of
        %   radial bins
        if ( iRadius + nRadiusTraveled(iVelo) ) > nRadius
            nRadiusDelta = nRadius - iRadius;
        else
            nRadiusDelta = nRadiusTraveled(iVelo);
        end

        %------------------------------------------------------------------
        % Get number of particles
        %------------------------------------------------------------------

        % Number of background particles in current bin
        nBG = bgMatrix(iVelo, iRadius);
        
        % Skip to next bin if the number of particles is below threshold
        if nBG < nMin
            continue
        end

        %------------------------------------------------------------------
        % Update background particle positions
        %------------------------------------------------------------------

        % Remove all background particles from current bin
        bgMatrix(iVelo, iRadius) = -nBG ...
            + bgMatrix(iVelo, iRadius);

        % Add all background particles to new radial bin
        bgMatrix(iVelo, iRadius + nRadiusDelta) = nBG ...
            + bgMatrix(iVelo, iRadius + nRadiusDelta);

        if debugBool
            % [DEBUG] Total loop counter
            loopCount = loopCount + 1;
        end
        
    end % Background velocity loop
    end % Background radial loop

    %----------------------------------------------------------------------
    % First pass to update velocities
    %----------------------------------------------------------------------
    
    for iRadius = (nRadius - 1) : -1 : 1
    %% Calculations per radial bin
    % Loop backwards to prevent counting particles twice

    for iVelo = nVelo : -1 : iFirstVelo
    %% Calculations per plasma velocity bin
    % Loop backwards as fast particles travel in front of slower particles

    %----------------------------------------------------------------------
    % Get number of particles
    %----------------------------------------------------------------------
    
    % Number of plasma particles in current bin
    nPlasma = plasmaMatrix(iVelo, iRadius);

    % Skip this velocity bin if the number of particles is
    %   lower than the threshold
    if nPlasma < nMin
        continue
    end
    
    %----------------------------------------------------------------------
    % Set number of traveled radial bins
    %----------------------------------------------------------------------

    % Restrict number of traveled bins to the total number of
    %   radial bins
    if ( iRadius + nRadiusTraveled(iVelo) ) > nRadius
        nRadiusDelta = nRadius - iRadius;
    else
        nRadiusDelta = nRadiusTraveled(iVelo);
    end
    
    % Total number of collided particles
    nCol = 0;

    % Loop through traversed bins
    for thisRadius = iRadius : iRadius + nRadiusDelta - 1
        %% Calculations per traversed radial bin

        %------------------------------------------------------------------
        % Get filled background particle bins
        %------------------------------------------------------------------

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

        % Skip to next radial bin if none of the velocity bin has enough
        %   particles
        if isempty(bgVeloArray)
            continue
        end
        
        % Update number of background particles array
        nBGVelo = nBGVelo(nBGVelo >= nMin);

        %------------------------------------------------------------------
        % Calculate velocity weight factors
        %------------------------------------------------------------------

        % Relative velocity weight factors
        veloWeightArray = (iVelo - bgVeloArray) ./ (iVelo + bgVeloArray);

        % Loop through the filled background velocities smaller than the
        %   plasma velocity
        for jVelo = 1 : numel(bgVeloArray)
        %% Calculations per bg velocity bin

            % Background velocity index
            iBGVelo = bgVeloArray(jVelo);

            %--------------------------------------------------------------
            % Collision rate per bin
            %--------------------------------------------------------------
            
            % Background density in current bin
            bgDensityBin = nBGVelo(jVelo) / binVolume(thisRadius);

            % Velocity weighted collision rate per bin
            colRateBin = bgDensityBin * sigma * radiusDelta ...
                         * veloWeightArray(jVelo);

            % Number of collided particles per bin
            nColBin = colRateBin * nPlasma;

            % Skip to next background velocity if not enough particles have
            %   collided
            if nColBin < nMin
                continue
            end
            
            % Normalization factor for average number of collisions per bin
            normK = 1;
            
            if nColBin > nBGVelo(jVelo)
                % Adjust the normalization factor if the number of collided
                %   particles is limited
                normK = nBGVelo(jVelo) / nColBin;
                
                % Limit the total number of collided particles to the
                %   number of background particles in this bin
                nColBin = nBGVelo(jVelo);
            end
            
%             % Update total number of collided particles
%             nCol = nCol + nColBin;
            
            %--------------------------------------------------------------
            % Calculate new velocities after collision
            %--------------------------------------------------------------

            % Calculate new plasma velocity after plasma-bg collision (eq. 6)
            %   * Only valid if veloPlasma >> veloBG
            iNewVelo = round( iVelo * (mass - massBG) ...
                                    / (mass + massBG) );

            % Calculate new bg velocity after plasma-bg collision (eq. 7)
            %   * Only valid if veloPlasma >> veloBG
            iNewBGVelo = round( 2 * mass * iVelo ...
                                 / (mass + massBG) );

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

            %--------------------------------------------------------------
            % Update background particle velocities
            %--------------------------------------------------------------

            % Remove collided bg particles from velo bin
            bgMatrix(iBGVelo, thisRadius) = -nColBin ...
                + bgMatrix(iBGVelo, thisRadius);

            % Add collided bg particles to new velo bin
            bgMatrix(iNewBGVelo, thisRadius) = nColBin ...
                + bgMatrix(iNewBGVelo, thisRadius);

            %--------------------------------------------------------------
            % Update plasma particle velocities and position
            %--------------------------------------------------------------

            % Remove collided plasma particles from the old bin
            plasmaMatrix(iVelo, iRadius) = -nColBin ...
                + plasmaMatrix(iVelo, iRadius);

            % Add collided plasma particles to the new bin
            plasmaMatrix(iNewVelo, iRadius + nRadiusDelta) = nColBin ...
                + plasmaMatrix(iNewVelo, iRadius + nRadiusDelta);

            % Increase number of collisions of bin by normalized collision
            %   rate
            plasmaMatrix(iNewVelo, iRadius + nRadiusDelta, kField) =  ...
                plasmaMatrix(iNewVelo, iRadius + nRadiusDelta, kField) ...
                + 1;

            if debugBool
                % [DEBUG] Total loop counter
                loopCount = loopCount + 1;
            end

        end % Background velocity loop

    end % Traversed radial bins loop
    
    %----------------------------------------------------------------------
    % Update uncollided plasma particle positions
    %----------------------------------------------------------------------
    
    % Number of uncollided plasma particles
    nPlasma = nPlasma - nCol;
    
    % Skip to next bin if the number of uncollided particles is below
    %   the threshold
    if nPlasma < nMin
        continue
    end
    
    % Prevent negative number of particles
    %   * Should not happen, built in for redundancy
    if nPlasma < 0
        nPlasma = 0;
    end

    % Remove all plasma particles from current bin
    plasmaMatrix(iVelo, iRadius) = -nPlasma ...
        + plasmaMatrix(iVelo, iRadius);

    % Add all plasma particles to new radial bin
    plasmaMatrix(iVelo, iRadius + nRadiusDelta) = nPlasma ...
        + plasmaMatrix(iVelo, iRadius + nRadiusDelta);

    end % Plasma velocity loop

    end % Radial loop

    %----------------------------------------------------------------------
    % Freeze particles in last radial bin
    %----------------------------------------------------------------------

    % Get all kinetic plasma particles in the last radial bin
    plasmaEnd = sum( plasmaMatrix(2:end, nRadius) );

    % If the number of kinetic pasma particles is above the threshold
    if plasmaEnd >= nMin
        % Remove plasma particles from the kinetic velocity bins
        plasmaMatrix(2:end, nRadius) = 0;

        % Add plasma particles to the static velocity bin
        plasmaMatrix(1, nRadius) = plasmaMatrix(1, nRadius) + plasmaEnd;
    end

    % Get all kinetic background particles in the last radial bin
    bgEnd = sum( bgMatrix(2:end, nRadius) );

    % If the number of kinetic bg particles is above the threshold
    if bgEnd >= nMin
        % Remove bg particles from the kinetic velocity bins
        bgMatrix(2:end, nRadius) = 0;

        % Add bg particles to the static velocity bin
        bgMatrix(1, nRadius) = bgMatrix(1, nRadius) + bgEnd;
    end
    
    

    %% Plot 1D propagation
    % Only for first angle (center of the plume)
    if iAngle == 1        
        % Only for specific times
        if (time(iTime) == 0.5E-6) || (time(iTime) == 1E-6) || ...
           (time(iTime) == 1.5E-6) || (round(time(iTime)*10^9) == 3000) || ...
           (time(iTime) == timeMax - timeDelta)

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
            
            % Increment plot index
            plotIndex = plotIndex + 1;
            
            % Plot 1D plasma propagation with separated collisions
            figure(figCollisions1D);
            subplot(1, 5, plotIndex);
            ylim([0 1E13]);
            title([num2str(time(iTime) * 1E6) ' \mus']);
            hold on;
            
            % Collision matrix
            colMatrix = round( plasmaMatrix(:, :, kField) );
            
            % Temporary matrix holding number of particles
            plasmaTemp = plasmaMatrix(:, :, nField);
            
            % Loop through number of collisions
            for iNCol = 0 : colMax
                
                % Loop through radial bins
                for iRadius = 1 : nRadius
                    % Get all bins with the current number of collisions
                    tempMatrix = plasmaTemp(colMatrix(:, iRadius) == iNCol, iRadius);

                    % Sum all particles with current number of collisions
                    tempRad(iRadius) = sum(tempMatrix);
                end

                % Plot the number of particles at current time step with
                %   current number of collisions
                plot(radius, tempRad);
            end
            
            % Plot legend only for the last plot
            if plotIndex == 1
                ylabel('Number of particles per radial bin');
            elseif plotIndex == 3
                xlabel('Target distance [m]');
            elseif plotIndex == 5
                legend('Uncollided', '1 collision', '2 collisions');
            end
            
        end
    end 

end % Temporal loop

end % Anglular loop

%% DEBUG

% Check final total amount of particles in matrices
endPlasma = sum(sum(plasmaMatrix(:,:,nField)));
endBG = sum(sum(bgMatrix));

% Check if the total number of plasma particles is conserved
if abs(startPlasma - endPlasma) < nMin
    disp('Total number of plasma particles was conserved.');
else
    disp(['Change in number of plasma particles detected: ' ...
            num2str(startPlasma - endPlasma, '%.3E')]);
%     error('Total number of plasma particles was NOT conserved.');
end

% Check if the total number of background particles is conserved
if abs(startBG - endBG) < nMin
    disp('Total number of background particles was conserved.');
else
    disp(['Change in number of background particles: ' ...
            num2str(startBG - endBG, '%.3E')]);
%     error('Total number of background particles was NOT conserved.');
end

% Temporary matrix holding number of particles
plasmaTemp = plasmaMatrix(:, :, nField);

% Total number of negative particles
nNeg = abs(sum(plasmaTemp(plasmaTemp < 0)));

% Check if there is a significant number of negative particles
if nNeg > nMin
    disp(['Significant number of negative plasma particles detected: ' ...
            num2str(nNeg, '%.3E')]);
end
    
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
