function particleMatrix = update2DMatrix( particleMatrix, nRadius, nVelo, nMin )
%UPDATEMATRIX Update a 2D particle matrix.
%   Note: updateMatrix is faster and should be used in most cases.

% Check is the input matrix is 2D
if ~ismatrix(particleMatrix)
    error('update2DMatrix can only be used for 2D matrices')
end

for iRadius = (nRadius - 1) : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

for iVelo = nVelo : -1 : 2
%% Calculations per plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Skip zero velocity

% Skip iteration if the number of particles is below the particle threshold
if particleMatrix(iVelo, iRadius) > nMin
    
    % Restrict number of traveled bins to the total number of radial bins,
    %   to prevent index out-of-bounds error
    if ( iRadius + iVelo - 1 ) > nRadius
        nRadiusDelta = nRadius - iRadius;
    else
        nRadiusDelta = iVelo - 1;
    end

    % Remove non-colliding particles to new position
    particleMatrix(iVelo, iRadius) = ...
        particleMatrix(iVelo, iRadius) - nParticle;

    % Add non-colliding particles to new position
    particleMatrix(iVelo, iRadius + nRadiusDelta) = ...
        particleMatrix(iVelo, iRadius + nRadiusDelta) + nParticle;

end

end % Velocity loop

end % Radius loop
    
end

