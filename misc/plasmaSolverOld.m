
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
