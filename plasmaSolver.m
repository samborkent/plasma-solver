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
%   velo    : Velocity
%   uc      : Unit cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializes

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
% Plot settings
%--------------------------------------------------------------------------

% Plot initial velocity distribution
plotVelDisInit      = true;

% Save the initial velocity distribution plot
saveVelDisInit      = false;

%--------------------------------------------------------------------------
% Computational restrictions
%--------------------------------------------------------------------------

% Minimal number of plasma particle per bin
minNPlasma = 1;

%--------------------------------------------------------------------------
% Model parameters
%--------------------------------------------------------------------------

% Fit parameter of the cosines power in the anular plasma particle
%   distibution function
cosPowerFit = 16;

%--------------------------------------------------------------------------
% Material settings
%--------------------------------------------------------------------------

% Unit cell
%   * See UnitCells
%   * Options: TiO2, SrTiO3, Li2O, LiO2, LiTi2O4, Li2TiO3, Li4Ti5O12
%   * Preface with UC.
uc = UC.TiO2;

%--------------------------------------------------------------------------
% Background gas parameters
%--------------------------------------------------------------------------

% Background gas pressure during depostion [Pa]
bgPressure      = 0.1E2;

% Temperature of background gas [K]
bgTemperature   = 300;      

%--------------------------------------------------------------------------
% Laser parameters
%--------------------------------------------------------------------------

% Laser spot width [m]
spotWidth       = 1E-3;

% Laser spot height [m]
spotHeight      = 2.1E-3;

% Ablation depth into the target for a single pulse [m]
ablationDepth   = 100E-9;

% Laser fluence [J / cm^2]
laserFluence    = 2.0;

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 9E-6;     % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.05;     % End position [m]
radiusDelta = 0.00005;  % Radial step size [m]

% Velocity limits
veloMax         = 2.5E4; % Maximal initial velocity
veloStepsize    = 5;     % Velocity step size

%% Calculations
% Everything below is calculated automatically based on user input above

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

%% Pre-allocations

% Plasma particle matrix
plasmaMatrix = zeros(nVelo, nRadius);

%% Angular plasma particle distribution

nParticleAngle = angularDistribution( angle, radiusDelta, ...
                                        cosPowerFit, nAtomAblated );

%% Initialize figures

% 1D propagation figure
figPropagation1D = figure('Name', 'Propagation 1D');
hold on;

%% Main program

sigma = pi * ( 2*PT.O.RADIUS + PT.Ti.RADIUS )^2;
mass = PT.Ti.MASS;
massBG = 2*PT.O.MASS;

for iAngle = 1 : 1 % (nAngle - 1)
    %% Calculations per angle
    
    % Skip angle if there are not enough particles present
    if nParticleAngle(iAngle) < minNPlasma
        continue
    end
    
    % Pre-allocate background gas particle matrix and fill based on
    %   calculated density
    bgMatrix = fillBGMatrix( bgDensity, nVelo, radius, angle, iAngle );
    
    % Calculate the initial particle velocity distribution
    nPlasmaVeloInit = initialVelocityDistribution( plotVelDisInit, ...
                        plotVelDisInit, velo, 1500, nUCAblated, uc, ...
                        1, energyLaser, 0, 0.6, ...
                        nParticleAngle(iAngle), 1 );
                    
    % Only plot and save initial particle velocity distribution for the
    %   center of the plume
    if iAngle == 1
        plotVelDisInit = false;
        saveVelDitInit = false;
    end

    % Fill into the plasma matrix
    plasmaMatrix(:, 1) = nPlasmaVeloInit(:, 1);

    for iTime = 1 : nTime
        %% Calculations per time step

        for iRadius = nRadius : -1 : 1
            %% Calculations per radial bin
            % Loop backwards to prevent counting particles twice

            for iVelo = nVelo : -1 : 1
                %% Calculations per plasma velocity bin
                % Loop backwards to prevent counting particles twice

                % Skip velocity bin if no radial bin is traveled within
                %   a time step
                if nRadiusTraveled(iVelo) == 0
                    continue
                end
                
                % Number of plasma particles in current bin
                nPlasma = plasmaMatrix(iVelo, iRadius);
                
                % Remove all particles from current bin
                plasmaMatrix(iVelo, iRadius) = 0;
                
                % Skip bin if the number of particles is lower than the
                %   minimum number of particles specified
                if nPlasma < minNPlasma                    
                    continue
                end

                % Restrict number of traveled bins to the total number of
                %   radial bins
                if ( iRadius + nRadiusTraveled(iVelo) ) > nRadius
                    nRadiusVelo = nRadius - iRadius;
                else
                    nRadiusVelo = nRadiusTraveled(iVelo);
                end
                
                %% TEST COLLISIONS
                
                % Collision probability of plasma particle with background
                %   gas particle
                colProb = timeDelta * velo(iVelo) * sigma * bgDensity;
                
                % If the collision probability is greater than 1, all
                %   particles collide
                if colProb > 1
                    colProb = 1;
                end
                
                colProb = colProb * 0.5;
                
                % Number of colliding particles
                nColPlasma = colProb * nPlasma;
                
                % Number of non-colliding particles
                nPlasma = nPlasma - nColPlasma;
                
                % Calculate new velocity after collision (eq. 6)
                iNewVelo = round( iVelo * ((mass - massBG) ...
                                        / (mass + massBG)) );
                
                % Prevent index out-of-range error
                if iNewVelo < 1
                    iNewVelo = 1;
                elseif iNewVelo > nVelo
                    iNewVelo = nVelo;
                end
                
                % Place colliding particles in new radial bin with altered
                %   velocity
                plasmaMatrix(iNewVelo, iRadius + nRadiusVelo) = nColPlasma ...
                    + plasmaMatrix(iNewVelo, iRadius + nRadiusVelo);
                
                % Place non-colliding particles in the new radial bin
                %   with unaltered velocity
                plasmaMatrix(iVelo, iRadius + nRadiusVelo) = nPlasma ...
                    + plasmaMatrix(iVelo, iRadius + nRadiusVelo);

            end % Velocity loop

        end % Radial loop
        
        %% Plot 1D propagation
        % Only for first angle (center of the plume)
        if iAngle == 1
            % Only for specific times
            if (time(iTime) == 0.1E-6) || (time(iTime) == timeMax - timeDelta) % || ...
%                (time(iTime) == 1E-6)   || (time(iTime) == 2E-6)   || ...
%                (time(iTime) == 3E-6)   || (time(iTime) == 6E-6)   || ...
%                (time(iTime) == timeMax - timeDelta)

                % Calculate total number of particles per radial bin
                nPlasmaRadius = nPlasmaParticlesPerRadius( radius, plasmaMatrix );

                % Get the envelope to smooth data
                [envNPlasmaRadius, ~] = envelope(nPlasmaRadius, 1, 'peak');

                % Normalization factor
                envNorm = nUCAblated / sum(envNPlasmaRadius);

                % Normalize envelope to conserve number of particles
                envNPlasmaRadius = envNPlasmaRadius .* envNorm;

                % Plot the 1D propagation
                figure(figPropagation1D);
                bar( radius, nPlasmaRadius, ...
                     'DisplayName', [num2str(time(iTime), 3) ' s'] );
%                 plot( radius, envNPlasmaRadius, ...
%                      'DisplayName', [num2str(time(iTime), 3) ' s'], ...
%                      'LineWidth', 2 );
            end
        end 

    end % Temporal loop

end % Anglular loop
    
%% Plot settings

% 1D propagation
figPropagation1D;
legend;
% ylim([0 2.5E14]);
title('1D propagation of Ti from TiO_2 target');
xlabel('Target distance [m]');
ylabel('Number of particles');

