function [particleMatrix, collisionMatrix, nLoops] = collisionCalculation( particleMatrix, ...
    nLoops, nMin, nMinBG, iRadiusRange, radiusDelta, nRadius, velo, veloDelta, nVelo, ...
    iVeloRange, iVeloZero, binVolume, colCS, mass, nSpecies )
%COLLISIONCALCULATION Perform all collision calculations.

% Initialize collision matrix
collisionMatrix = zeros(nSpecies, 1, nVelo, nRadius);

%-------------------------------------------------------------------------------
% Pre-calculations
%-------------------------------------------------------------------------------
                             
% Calculate collision rate matrix
%   * Density of background particles, times length of one radial bin,
%       time the collision cross-section
%   * Prevents having to calculate the collision rate millions of times,
%       instead calculate the collision rate for every position at once
%       each time step
collisionRate = ( 1 ./ reshape(binVolume, 1, 1, 1, nRadius) ) ...
                .* radiusDelta .* colCS;
            
% Matrix containing all possible velocity weights for collision calculation
%   * Prevents having to calculate many times, instead just look up value
%       from matrix
veloWeight = abs( (velo' - velo) ./ (abs(velo') + abs(velo)) );

for iRadius = flip( iRadiusRange( any( particleMatrix(2:end, :, :, 1:end-1) > nMin, [1 2 3] ) ) )
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice (flip)
% Skip last bin as particles cannot move further (1:end-1)
% Only loop through radial bins where any of the bins inside have a number
%   of plasma particles above the threshold

for iVelo = flip( iVeloRange( any( particleMatrix(2:end, :, :, iRadius) > nMin, [1 2] ) ) )
%% Calculations per non-oxidizing plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles (flip)
% Only loop through velocity bins where any of the bins inside have a number
%   of plasma particles above the threshold

% Skip zero velocity as no distance is traversed
if iVelo == iVeloZero
    continue
end

% Save number of plasma particles before any collisions occured
nPlasmaTemp = particleMatrix(2:end, :, iVelo, iRadius);

%-------------------------------------------------------------------------------
% Determine number of traveled bins without collisions
%-------------------------------------------------------------------------------

% If velocity if positive
if iVelo > iVeloZero
    % If the number of forwards traveled radial bins exceeds the last bin
    if ( iRadius + (iVelo - iVeloZero) ) > nRadius
        % Limit to distance to last bin
        nRadiusDelta = nRadius - iRadius;
    else
        % Set number of traveled bins based on velocity
        nRadiusDelta = iVelo - iVeloZero;
    end
% If velocity if negative
elseif iVelo < iVeloZero
    % If the number of backwards traveled radial bins exceeds the first bin
    if  ( iRadius + (iVelo - iVeloZero) ) < 1
        % Limit to distance to first bin
        nRadiusDelta = iRadius - 1;
    else
        nRadiusDelta = abs(iVelo - iVeloZero);
    end
end

% If no distance is traveled skip to next bin
if nRadiusDelta == 0
    continue
end

for jRadius = 0 : (nRadiusDelta - 1)
%% Calculations per traversed radial bin
% Loop from current bin to one bin before final bin
% Collisions occur first in the bin closest to starting bin

% Break loop if no more particles are left to collide
if all( nPlasmaTemp < nMin, 'all' )
    break
end

% Index of the current radius: starting bin + traveled bin
if iVelo > iVeloZero
    thisRadius = iRadius + jRadius;
elseif iVelo < iVeloZero
    thisRadius = iRadius - jRadius;
end

% Find the filled background particle velocity bins at this radius
if iVelo > iVeloZero
    iVeloBGArray = iVeloRange(1:iVelo-1);
    iVeloBGArray = iVeloBGArray( any(particleMatrix(1, :, 1:iVelo-1, thisRadius) > nMinBG, 2) );
elseif iVelo < iVeloZero
    iVeloBGArray = iVeloRange(iVelo+1:end);
    iVeloBGArray = iVeloBGArray( any(particleMatrix(1, :, iVelo+1:end, thisRadius) > nMinBG, 2) );
end

for iVeloBG = iVeloBGArray
%% Calculations per background velocity bin
% Only include background velocities smaller than the plasma velocity
% From slowest to highest as collisions with slow moving background
%   particles are more probable

% Break loop if no more particles are left to collide
if all( nPlasmaTemp < nMin, 'all' )
    break
end

% Number of background particles
nBG = sum(particleMatrix(1, :, iVeloBG, thisRadius), 2);

% If number of background particles is below threshold
if nBG < nMinBG
    continue
end

% Number of collisions array
nCol = collisionRate(:, 1, 1, thisRadius) ...   % Pre-calculated collision rate
       .* nBG...                                % Number of background particles
       .* veloWeight(iVelo, iVeloBG) ...        % Difference in velocity weight
       .* nPlasmaTemp;                          % Number of plasma particles left

% Skip iteration if the number of collisions is smaller than the threshold
if all( nCol < nMinBG, 'all' )
    continue
end

% If the number of collisions is greater than the number of background particles
if sum( nCol, 'all' ) > nBG
    % Normalize number of collisions to number of background particles
    %   * Minus 1 to avoid faulty comparisons
    nCol = (nCol ./ sum(nCol, 'all')) .* (nBG - 1);
end

% Update number of non-collided plasma particles
nPlasmaTemp = nPlasmaTemp - nCol;

%-------------------------------------------------------------------------------
% Calculate new velocity after collision
%-------------------------------------------------------------------------------

% Calculate the new velocity of plasma particles after collision
iNewVelo = iVeloZero + round( ( velo(iVelo).*(mass(2:end) - mass(1)) ...
                                + 2*mass(1)*velo(iVeloBG) ) ...
                              ./ ( mass(2:end) + mass(1) ) ./ veloDelta );

% Calculate the new velocity background gas particles after collision
iNewVeloBG = iVeloZero + round( ( velo(iVeloBG).*(mass(1) - mass(2:end)) ...
                                  + 2*mass(2:end)*velo(iVelo) ) ...
                                ./ ( mass(1) + mass(2:end) ) ./ veloDelta );

% Prevent index out-of-range error
iNewVelo(iNewVelo < 1) = 1;
iNewVelo(iNewVelo > nVelo) = nVelo;
iNewVeloBG(iNewVeloBG < 1) = 1;
iNewVeloBG(iNewVeloBG > nVelo) = nVelo;

%-------------------------------------------------------------------------------
% Calculate new radial positions after collision
%-------------------------------------------------------------------------------

% New radial index of collided plasma particles
nNewRadiusDelta = round( velo(iNewVelo) ...
                         .* (nRadiusDelta - jRadius) ...
                         ./ (nRadiusDelta * veloDelta) );
                       
% New radial index of collided background particles
nNewRadiusDeltaBG = round( velo(iNewVeloBG) ...
                           .* (nRadiusDelta - jRadius) ...
                           ./ (nRadiusDelta * veloDelta) );

% Prevent index out-of-range error
nNewRadiusDelta(thisRadius + nNewRadiusDelta < 1) = ...
    1 - thisRadius;
nNewRadiusDeltaBG(thisRadius + nNewRadiusDeltaBG < 1) = ...
    1 - thisRadius;
nNewRadiusDelta(thisRadius + nNewRadiusDelta > nRadius) = ...
    nRadius - thisRadius;
nNewRadiusDeltaBG(thisRadius + nNewRadiusDeltaBG > nRadius) = ...
    nRadius - thisRadius;

%-------------------------------------------------------------------------------
% Update collided particles
%-------------------------------------------------------------------------------

% % Remove plasma particles from initial position before collision
particleMatrix(2:end, :, iVelo, iRadius) = ...
    particleMatrix(2:end, :, iVelo, iRadius) - nCol;

% Remove background particles from initial position before collision
%   * The total number of collided particles get removed from the
%       non-collided and collided background gas particles bins based on
%       the ratio between non-collided and collided particles.
particleMatrix(1, :, iVeloBG, thisRadius) = ...
    particleMatrix(1, :, iVeloBG, thisRadius) ...
    - ( particleMatrix(1, :, iVeloBG, thisRadius) ...
        ./ sum(particleMatrix(1, :, iVeloBG, thisRadius), 2) ) ...
      * sum(nCol, 'all');

% Sum number of collisions over number of collisions per particle
nColSum = sum(nCol, 2);

% Loop through plasma species
for iSpecies = 2 : nSpecies
    % Add plasma particles to new position after collision
    collisionMatrix(iSpecies, 1, iNewVelo(iSpecies-1), thisRadius + nNewRadiusDelta(iSpecies-1)) = ...
        collisionMatrix(iSpecies, 1, iNewVelo(iSpecies-1), thisRadius + nNewRadiusDelta(iSpecies-1)) ...
        + nColSum(iSpecies-1);

    % Add background particles to new position after collision
    collisionMatrix(1, 1, iNewVeloBG(iSpecies-1), thisRadius + nNewRadiusDeltaBG(iSpecies-1)) = ...
        collisionMatrix(1, 1, iNewVeloBG(iSpecies-1), thisRadius + nNewRadiusDeltaBG(iSpecies-1)) ...
        + nColSum(iSpecies-1);
end

% Increment number of total loops
nLoops = nLoops + 1;

end % Background velocity loop

end % Traveled distance loop

end % Velocity loop

end % Radius loop

end

