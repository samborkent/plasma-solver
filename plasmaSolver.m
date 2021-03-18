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

% Minimal number of particle per bin
nMin = 1;

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
timeMax     = 7E-6;     % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 5E-5;  % Radial step size [m]

% Velocity limits
veloMax         = 2.5E4; % Maximal initial velocity
veloStepsize    = 5;    % Velocity step size

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

nParticleAngle = angularDistribution( angle,        ...
                                      radiusDelta,  ...
                                      cosPowerFit,  ...
                                      nAtomAblated );

%% Initialize figures

% 1D propagation figure
figPropagation1D = figure('Name', 'Propagation 1D');
hold on;

%% Main program

sigma = pi * ( 2*PT.O.RADIUS + PT.Ti.RADIUS )^2;
sigmaBG = pi * ( 4*PT.O.RADIUS )^2;
mass = PT.Ti.MASS;
massBG = 2*PT.O.MASS;

for iAngle = 1 : 1 % (nAngle - 1)
    %% Calculations per angle
    
    % Skip angle if there are not enough particles present
    if nParticleAngle(iAngle) < nMin
        continue
    end
    
    % Pre-allocate background gas particle matrix and fill based on
    %   calculated density
    [bgMatrix, binVolume] = fillBGMatrix( bgDensity,    ...
                                          nVelo,        ...
                                          radius,       ...
                                          angle,        ...
                                          iAngle );
    
    % Calculate the initial particle velocity distribution
    nPlasmaVeloInit = initialVelocityDistribution( plotVelDisInit,          ...
                                                   saveVelDisInit,          ...
                                                   velo,                    ...
                                                   1500,                    ...
                                                   nUCAblated,              ...
                                                   uc,                      ...
                                                   1,                       ...
                                                   energyLaser,             ...
                                                   0,                       ...
                                                   0.6,                     ...
                                                   nParticleAngle(iAngle),  ...
                                                   1 );
                    
    % Only plot and save initial particle velocity distribution for the
    %   center of the plume
    if iAngle == 1
        plotVelDisInit = false;
        saveVelDitInit = false;
    end

    % Fill into the plasma matrix
    plasmaMatrix(:, 1) = nPlasmaVeloInit(:, 1);
    
    % Check initial total amount particles in matrix [DEBUG]
    startPlasma = sum(sum(plasmaMatrix));
    startBG = sum(sum(bgMatrix));

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
                
                % Skip bin if the number of particles is lower than the
                %   minimum number of particles specified
                if nPlasma < nMin               
                    continue
                end
                
                % Remove all plasma particles from current bin
                plasmaMatrix(iVelo, iRadius) = 0;
                
                % Restrict number of traveled bins to the total number of
                %   radial bins
                if ( iRadius + nRadiusTraveled(iVelo) ) > nRadius
                    nRadiusDelta = nRadius - iRadius;
                else
                    nRadiusDelta = nRadiusTraveled(iVelo);
                end
                
                % Number of kinetic background particles in current bin
                nBGKinetic = bgMatrix(iVelo, iRadius);
                
                % Remove all kinetic particles from current bin
                bgMatrix(iVelo, iRadius) = 0;
                
                % Mean background density of traversed path
                meanBGDensity = ...
                    sum(sum( bgMatrix(:, iRadius+1:iRadius+nRadiusDelta) )) ...
                      / sum( binVolume(iRadius+1:iRadius+nRadiusDelta) );
                
                % Mean collision probability
                colProb = meanBGDensity * sigma * radiusDelta;
                
                % Mean collision probability per traversed bin
                colProbPerBin = colProb / nRadiusDelta;
                
                % Reset total number of collided particles [DEBUG]
                nColPlasma = 0;

                % Loop through traversed bins
                for iRadiusDelta = 1 : nRadiusDelta
                    
                    % Current radial bin index
                    thisRadius = iRadius + iRadiusDelta;
                    
                    % Find velocity bins containing particles below the
                    %   plasma velocity
                    bgVeloArray = find( ...
                        bgMatrix(1:iVelo-1, thisRadius) );
                    
                    % Loop through background velocities smaller than
                    %   plasma velocity
                    for jVelo = 1 : numel(bgVeloArray)
                        
                        % Background velocity index
                        iBGVelo = bgVeloArray(jVelo);
                        
                        % Skip to next velocity if number of background
                        %   particles is smaller than minimum amount
                        if bgMatrix(iBGVelo, thisRadius) < nMin
                            continue
                        end
                        
                        % Velocity weighted collision probability per bin
                        colProbPerBinNorm = colProbPerBin ...
                                            * (iVelo - iBGVelo) ...
                                            / (iVelo + iBGVelo);
                        
                        % Number of collided particles per bin
                        nColPlasmaPerBin = colProbPerBinNorm * nPlasma;
                        
                        % If enough particles have collided
                        if nColPlasmaPerBin >= nMin
                            
                            % Calculate new plasma velocity after collision (eq. 6)
                            iNewVelo = round( iVelo * ((mass - massBG) ...
                                                    /  (mass + massBG)) );

                            % Calculate new background velocity after plasma-bg collision (eq. 7)
                            iNewBGVelo = round( 2 * mass * iVelo / (mass + massBG) );

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
                            
                            % Update total number of collided particles [DEBUG]
                            nColPlasma = nColPlasma + nColPlasmaPerBin;
                            
                            % Update number of non-colliding particles
                            nPlasma = nPlasma - nColPlasmaPerBin;
                            
                            % Remove collided bg particles from velo bin
                            bgMatrix(iBGVelo, thisRadius) = ...
                                bgMatrix(iBGVelo, thisRadius) ...
                                - nColPlasmaPerBin;
                            
                            % Add collided bg particles to new velo bin
                            bgMatrix(iNewBGVelo, thisRadius) = ...
                                bgMatrix(iNewBGVelo, thisRadius) ...
                                + nColPlasmaPerBin;
                            
                            % Add collided plasma particle to new velo bin
                            plasmaMatrix(iNewVelo, thisRadius) = ...
                                plasmaMatrix(iNewVelo, thisRadius) ...
                                + nColPlasmaPerBin;
                            
                        end
                        
                        % If the number of collided plasma particles is
                        %   less than the number of background particles,
                        %   all plasma particles will be scattered in this
                        %   bin, so skip the rest of the loop
                        if nColPlasmaPerBin < bgMatrix(iBGVelo, thisRadius)
                            break
                        end
                        
                    end % Background velocity loop
                    
                end % Traversed radial bins loop
                
                % Place non-colliding particles in new radial bin
                %   with unaltered velocity
                plasmaMatrix(iVelo, iRadius + nRadiusDelta) = nPlasma ...
                    + plasmaMatrix(iVelo, iRadius + nRadiusDelta);
                
                % Place kinetic background particles in new radial bin
                bgMatrix(iVelo, iRadius + nRadiusDelta) = nBGKinetic ...
                    + bgMatrix(iVelo, iRadius + nRadiusDelta);
                
%                 % Number of static background particles in radial bin
%                 nBGStatic = bgMatrix(1, iRadius);
%                 
%                 % Initialize number of kinetic background particles
%                 nBGKinetic = 0;
%                 
%                 % If velocity is non-zero
%                 if iVelo > 1
%                     % Number of background particles in velocity bin
%                     nBGKinetic = bgMatrix(iVelo, iRadius);
%                     
%                     % Remove all background particles from velocity bin
%                     bgMatrix(iVelo, iRadius) = 0;
%                 end
%                 
%                 %% TEST COLLISIONS plasma - bg
%                 
%                 % Calculate new plasma velocity after collision (eq. 6)
%                 iNewVelo = round( iVelo * ((mass - massBG) ...
%                                         / (mass + massBG)) );
%                                     
%                 % Calculate new background velocity after plasma-bg collision (eq. 7)
%                 jNewVelo = round( 2 * mass * iVelo / (mass + massBG) );
%                 
%                 % Prevent index out-of-range error
%                 if iNewVelo < 1
%                     iNewVelo = 1;
%                 elseif iNewVelo > nVelo
%                     iNewVelo = nVelo;
%                 end
%                 if jNewVelo < 1
%                     jNewVelo = 1;
%                 elseif jNewVelo > nVelo
%                     jNewVelo = nVelo;
%                 end
%                 
%                 % Collision probability of plasma particle with static
%                 %   background gas particles
%                 colProbPlasma = timeDelta * velo(iVelo) * sigma ...
%                                 * (nBGStatic / binVolume(iRadius));
%                 
%                 % If the collision probability is greater than 1, all
%                 %   particles collide
%                 if colProbPlasma > 1
%                     colProbPlasma = 1;
%                 end
%                 
%                 % Number of colliding particles
%                 nColPlasma = colProbPlasma * nPlasma;
%                 
%                 % If enough particles have collided
%                 if nColPlasma >= nMin
%                     % Number of non-colliding particles
%                     nPlasma = nPlasma - nColPlasma;
%                     
%                     % Place colliding particles in new radial bin with altered
%                     %   velocity
%                     plasmaMatrix(iNewVelo, iRadius + nRadiusDelta) = nColPlasma ...
%                         + plasmaMatrix(iNewVelo, iRadius + nRadiusDelta);
%                     
%                     % Remove collided particles from current bg bin
%                     bgMatrix(1, iRadius) = -nColPlasma ...
%                         + bgMatrix(iVelo, iRadius);
%                     
%                     % Add them to new bg velocity bin
%                     bgMatrix(jNewVelo, iRadius) = nColPlasma ...
%                         + bgMatrix(jNewVelo, iRadius);
%                 end
%                 
%                 % If enough background particles are in the velocity bin
%                 if nBGKinetic >= nMin
%                     % Place kinetic background particles in new radial bin
%                     bgMatrix(iVelo, iRadius + nRadiusDelta) = nBGKinetic ...
%                         + bgMatrix(iVelo, iRadius + nRadiusDelta);
%                 end
%                 
%                 if colProbPlasma > 0
%                     disp( ['iTime: ' num2str(iTime)] );
%                     disp( ['iRadius: ' num2str(iRadius)] );
%                     disp( ['iVelo: ' num2str(iVelo)] );
%                     disp( ['nRadiusDelta: ' num2str(nRadiusDelta)] );
%                     disp( ['colProbPlasma: ' num2str(colProbPlasma)] );
%                     disp( ['nBG: ' num2str(nBGStatic)] );
%                     disp( ['nColPlasma: ' num2str(nColPlasma)] );
%                     disp( ['nPlasma: ' num2str(nPlasma)] );
%                     disp( ['iNewVelo: ' num2str(iNewVelo)] );
%                     disp( ['jNewVelo: ' num2str(jNewVelo)] );
%                     disp( ' ' );
%                 end

            end % Plasma velocity loop

        end % Radial loop
        
        %% Plot 1D propagation
        % Only for first angle (center of the plume)
        if iAngle == 1
            % Only for specific times
            if (time(iTime) == timeDelta) || (time(iTime) == 0.5E-6) || ...
               (time(iTime) == 1E-6)      || (time(iTime) == 2E-6)   || ...
               (time(iTime) == 3E-6)      || (time(iTime) == 6E-6)   || ...
               (time(iTime) == timeMax - timeDelta)

                % Calculate total number of particles per radial bin
                nPlasmaRadius = nPlasmaParticlesPerRadius( radius, plasmaMatrix );
                
                % Smooth the number of particle data
%                 nPlasmaRadius = smooth(nPlasmaRadius, veloStepsize);

                % Normalization factor
                plasmaNorm = nUCAblated / sum(nPlasmaRadius);

                % Normalize envelope to conserve number of particles
                nPlasmaRadius = nPlasmaRadius .* plasmaNorm;

                % Plot the 1D propagation
                figure(figPropagation1D);
                bar( radius, nPlasmaRadius, ...
                     'DisplayName', [num2str(time(iTime), 3) ' s'], ...
                     'LineStyle', 'none');
%                 plot( radius, nPlasmaRadius, ...
%                      'DisplayName', [num2str(time(iTime), 3) ' s'], ...
%                      'LineWidth', 2 );
            end
        end 

    end % Temporal loop

end % Anglular loop

%% DEBUG

% Check final total amount of particles in matrices [DEBUG]
endPlasma = sum(sum(plasmaMatrix));
endBG = sum(sum(bgMatrix));

% Check if total number of plasma particles is conserved
if abs(startPlasma - endPlasma) < nMin
    disp('Total number of plasma particles was conserved.');
else
    disp(['Change in number of plasma particles: ' ...
            num2str(startPlasma - endPlasma)]);
    error('Total number of plasma particles was NOT conserved.');
end

% Check if total number of background particles is conserved
if abs(startBG - endBG) < nMin
    disp('Total number of background particles was conserved.');
else
    disp(['Change in number of background particles: ' ...
            num2str(startBG - endBG)]);
    error('Total number of background particles was NOT conserved.');
end
    
%% Plot settings

% 1D propagation
figPropagation1D;
legend;
xlim([0 0.06]);
% ylim([0 2.5E14]);
title('1D propagation of Ti from TiO_2 target');
xlabel('Target distance [m]');
ylabel('Number of particles');

