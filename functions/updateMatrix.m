function particleMatrix = updateMatrix( particleMatrix, ...
                                        nRadius, nVelo, nMin )
%UPDATEMATRIX Update non-colliding particle matrix

% Number of species in particle matrix
nSpecies = size(particleMatrix, 2);

% Get rid of dimensions of size 1
if size(particleMatrix, 2) == 1
    particleMatrix = squeeze(particleMatrix);
end

%% If the matrix has 2 dimensions
if ismatrix(particleMatrix)

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
    nParticle = particleMatrix(iVelo, iRadius);

    % Skip iteration if the number of particles is below the threshold
    if nParticle < nMin
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
    particleMatrix(iVelo, iRadius) = ...
        particleMatrix(iVelo, iRadius) - nParticle;

    % Add non-colliding particles to new position
    particleMatrix(iVelo, iRadius + nRadiusDelta) = ...
        particleMatrix(iVelo, iRadius + nRadiusDelta) + nParticle;

    end % Velocity loop

    end % Radius loop
    
%% If the matrix has more than 2 dimensions
else

    for iRadius = (nRadius - 1) : -1 : 1
    %% Calculations per radial bin
    % Loop backwards to prevent counting particles twice
    % Skip last bin as particles cannot move further
    
    % Calculate sum of all particles in radial bin
    nParticleRadius = sum( particleMatrix(:, :, iRadius) );
    
    % If no particles are present in radial bin, skip to next radial bin
    if sum(nParticleRadius) < (nMin * nSpecies)
        continue
    else
        % Index of non-empty species
        iSpeciesRadius = find(nParticleRadius >= nMin);
    end

    for iVelo = nVelo : -1 : 2
    %% Calculations per plasma velocity bin
    % Loop backwards as fast particles travel in front of slower particles
    % Only loop through velocities that traverse at least one radial bin per
    %   time step
    
    %--------------------------------------------------------------------------
    % Get plasma particle variables
    %--------------------------------------------------------------------------
    
    % Number of plasma particles in current bin
    nParticle = particleMatrix(iVelo, iSpeciesRadius, iRadius);
    
    % Logical array holding locations of velocity bins with enough particles
    iParticle = nParticle >= nMin;
    
    % Remove empty plasma species
    nParticle = nParticle(iParticle);
    
    % If no filled bins remain, skip to next radial bin
    if isempty(nParticle)
        continue
    end
    
    % Indices of non-empty plasma species bins
    iSpeciesVelo = iSpeciesRadius(iParticle);

    %--------------------------------------------------------------------------
    % Set number of traveled radial bins
    %--------------------------------------------------------------------------

    % Restrict number of traveled bins to the total number of
    %   radial bins
    if ( iRadius + iVelo - 1 ) > nRadius
        nRadiusDelta = nRadius - iRadius;
    else
        nRadiusDelta = iVelo - 1;
    end

    %--------------------------------------------------------------------------
    % Update non-colliding particles
    %--------------------------------------------------------------------------
    
    % Remove non-colliding particles to new position
    particleMatrix(iVelo, iSpeciesVelo, iRadius) = ...
        particleMatrix(iVelo, iSpeciesVelo, iRadius) - nParticle;

    % Add non-colliding particles to new position
    particleMatrix(iVelo, iSpeciesVelo, iRadius + nRadiusDelta) = ...
        particleMatrix(iVelo, iSpeciesVelo, iRadius + nRadiusDelta) + nParticle;

    end % Velocity loop

    end % Radius loop
    
end % If ismatrix

end

