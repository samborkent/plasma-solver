function nParticles = nParticlesPerRadius( radius, particleMatrix )
%NPLASMARADIUS Calculate total number of plasma particles per radial bin

% Number of radial bins
nRadius = numel(radius);

% Pre-allocate
nParticles = zeros(1, nRadius);

% Loop through radial bins
for iRadius = 1 : nRadius
    % Sum particles in all velocity bins within current radial bin
    nParticles(iRadius) = sum( particleMatrix(:, iRadius) );
end
    
end

