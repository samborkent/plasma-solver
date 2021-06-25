function particleMatrix = update2DMatrixNegVelo( particleMatrix, ...
                                                 nRadius, nVelo, iVeloZero, nMin )
%UPDATEMATRIX Update a 2D particle matrix with support for negative velocities.
%   Note: updateMatrix is faster and should be used in most cases.

% Get dimensions of particle matrix
dims = size( particleMatrix );

% Remove dimentions with size 1
particleMatrix = squeeze( particleMatrix );

% Check is the input matrix is 2D
if ~ismatrix(particleMatrix)
    error('update2DMatrixNegVelo can only be used for 2D matrices')
end

% Matrix with new positions
newMatrix = particleMatrix.*0;

for iRadius = nRadius : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

for iVelo = nVelo : -1 : 1
%% Calculations per plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles

% Skip zero velocity
if iVelo == iVeloZero
    continue
end

% Skip iteration if the number of particles is below the particle threshold
if particleMatrix(iVelo, iRadius) > nMin
    
    % Number of particles
    nParticle = particleMatrix(iVelo, iRadius);
    
    % Restrict number of traveled bins to the total number of radial bins,
    %   to prevent index out-of-bounds error
    if iVelo > iVeloZero
        if iRadius + iVelo > nRadius
            nRadiusDelta = nRadius - iRadius;
        else
            nRadiusDelta = iVelo - iVeloZero;
        end
    elseif iVelo < iVeloZero
        if iRadius + (iVelo - iVeloZero) <= 0
            nRadiusDelta = 1 - iRadius;
        else
            nRadiusDelta = iVelo - iVeloZero;
        end
    end
    
    % Remove non-colliding particles from old position
    particleMatrix(iVelo, iRadius) = 0;

    % Add non-colliding particles to new position
    newMatrix(iVelo, iRadius + nRadiusDelta) = ...
        newMatrix(iVelo, iRadius + nRadiusDelta) + nParticle;
end

end % Velocity loop

end % Radius loop

% Add particles to new position
particleMatrix = particleMatrix + newMatrix;

% Reshape matrix to original dimensions
particleMatrix = reshape( particleMatrix, dims );
    
end

