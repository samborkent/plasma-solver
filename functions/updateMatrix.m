function updatedMatrix = updateMatrix( particleMatrix, nRadius, nVelo, ...
                                       iFirstVelo, nMin, nRadiusTraveled )
%UPDATEMATRIX Update non-colliding particle matrix

% Pre-allocate updated matrix
updatedMatrix = particleMatrix.*0;

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
% Get number of particles
%--------------------------------------------------------------------------

% Number of plasma particles in current bin
nParticle = particleMatrix(iVelo, iRadius);

% Skip iteration if the number of particles is below the threshold
if nParticle < nMin
    continue
end

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

%--------------------------------------------------------------------------
% Update non-colliding particles
%--------------------------------------------------------------------------

% Remove non-colliding particles to new position
updatedMatrix(iVelo, iRadius) = ...
    updatedMatrix(iVelo, iRadius) - nParticle;

% Add non-colliding particles to new position
updatedMatrix(iVelo, iRadius + nRadiusDelta) = ...
    updatedMatrix(iVelo, iRadius + nRadiusDelta) + nParticle;

end % Velocity loop

end % Radius loop

% Update matrix
updatedMatrix = particleMatrix + updatedMatrix;

end

