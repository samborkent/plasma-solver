function colMatrix = updateCol( colMatrix, nRadius, nVelo, nMin )
%UPDATECOl Update collision matrix

for iRadius = (nRadius - 1) : -1 : 1
%% Calculations per radial bin
% Loop backwards to prevent counting particles twice
% Skip last bin as particles cannot move further

for iVelo = nVelo : -1 : 2
%% Calculations per plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Only loop through velocities that traverse at least one radial bin
%   per time step

%----------------------------------------------------------------------
% Get number of particles
%----------------------------------------------------------------------

% Number of plasma particles in current bin
nParticle = colMatrix(:, iVelo, iRadius);

if sum(nParticle) < nMin
    continue
end

%----------------------------------------------------------------------
% Set number of traveled radial bins
%----------------------------------------------------------------------

% Restrict number of traveled bins to the total number of
%   radial bins
if ( iRadius + iVelo - 1 ) > nRadius
    nRadiusDelta = nRadius - iRadius;
else
    nRadiusDelta = iVelo - 1;
end

%----------------------------------------------------------------------
% Update non-colliding particles
%----------------------------------------------------------------------

% Remove non-colliding particles to new position
colMatrix(:, iVelo, iRadius) = colMatrix(:, iVelo, iRadius) - nParticle;

% Add non-colliding particles to new position
colMatrix(:, iVelo, iRadius + nRadiusDelta) = ...
    colMatrix(:, iVelo, iRadius + nRadiusDelta) + nParticle;

end % Velocity loop

end % Radius loop

end

