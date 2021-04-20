function nParticles = nParticlesPerRadius( radius, particleMatrix )
%NPLASMARADIUS Calculate total number of plasma particles per radial bin

% Remove dimensions with size 1
squeezeMatrix = squeeze( particleMatrix );

% Number of radial bins
nRadius = numel(radius);

% Pre-allocate
nParticles = zeros(1, nRadius);

% Loop through radial bins
for iRadius = 1 : nRadius
    % Sum particles in all velocity bins within current radial bin
    nParticles(iRadius) = sum( squeezeMatrix(:, iRadius) );
end
    
end

