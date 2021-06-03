function nParticleAngle = angularDistribution( angle, angleDelta, radiusDelta, ...
                                               cosPowerFit, nAtomAblated )
%ANGULARDISTRIBUTION Calculate the angular plasma particle distribution

% Number of angles
nAngle = numel(angle);

% Pre-allocate memory
nParticleAngle = zeros(1, nAngle);

% % Compute initial angular plasma particle distribution (eq. 5)
% for iAngle = 1 : nAngle
%     nParticleAngle(iAngle) = cosd(angle(iAngle)).^cosPowerFit;
% end

% Compute initial angular plasma particle distribution (eq. 5)
for iAngle = 1 : (nAngle - 1)
    nParticleAngle(iAngle) = ((4/3)*pi * radiusDelta^3) ...
        .* ( -cosd(angle(iAngle + 1)) + cosd(angle(iAngle)) ) ...
        .* cosd(angle(iAngle)).^cosPowerFit;
end

% Normalize and multiply by the total number of ablated atoms
% nParticleAngle = ( nParticleAngle ./ sum(nParticleAngle) ) ...
%                  .* ( nAtomAblated * deg2rad(angleDelta) );
nParticleAngle = ( nParticleAngle ./ sum(nParticleAngle) ) ...
                 .* nAtomAblated;

end

