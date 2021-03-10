function nPlasmaRadius = nPlasmaParticlesPerRadius( radius, plasmaMatrix )
%NPLASMARADIUS Calculate total number of plasma particles per radial bin

% Number of radial bins
nRadius = numel(radius);

% Pre-allocate
nPlasmaRadius = zeros(1, nRadius);

% Loop through radial bins
for iRadius = 1 : nRadius
    % Sum particles in all velocity bins within current radial bin
    nPlasmaRadius(iRadius) = sum( plasmaMatrix(:, iRadius) );
end
    
end

