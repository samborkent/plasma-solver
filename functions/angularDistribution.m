function nParticleAngle = angularDistribution( angle, radiusDelta, ...
    cosPowerFit, nAtomAblated )
%ANGULARDISTRIBUTION Calculate the angular plasma particle distribution

nAngle = numel(angle);

% Pre-allocate memory
nParticleAngle = zeros(1, nAngle - 1);

% Compute initial angular plasma particle distribution (eq. 5)
for iAngle = 1 : (nAngle - 1)
    nParticleAngle(iAngle) = ((4/3)*pi * radiusDelta^3) ...
        .* abs( cosd(angle(iAngle + 1)) - cosd(angle(iAngle)) ) ...
        .* cosd(angle(iAngle)).^cosPowerFit;
end

% Normalize and multiply by the total number of ablated atoms
nParticleAngle = (nParticleAngle .* nAtomAblated) ./ sum(nParticleAngle);

end

