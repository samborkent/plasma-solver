function colProb = collisionProbability( nParticle, ...
                                         startRadius, ...
                                         iVelo, ...
                                         radiusDelta, ...
                                         nRadiusTraveled, ...
                                         scatterCS, ...
                                         binVolumeArray, ...
                                         collisionMatrix )
%COLLISIONPROBABILITY Calculate the collision probability

%% Pre-allocate arrays
% Array of collision probability per radial bin traveled
colProb = zeros(1, nRadiusTraveled);

% Array of maximum number of collided particles
maxNColArray = ones(1, nRadiusTraveled) .* nParticle;

%% Main function

% Loop through the radius bins traveled by the propagating species
for iRadiusTraveled = 1 : nRadiusTraveled
    
    % Maximum number of collisions given by the sum of collision medium
    %   particles over all velocity bins in current radial bin
    maxNCol = sum(collisionMatrix(:, startRadius + iRadiusTraveled));
    
    % Add to array holding the maximum number of collisions per radial bin
    if nParticle > sum(collisionMatrix(:, startRadius + iRadiusTraveled))
        maxNColArray(iRadiusTraveled) = maxNCol;
    end
    
    % Find indices of non-empty velocity bins of the collision medium
    jVeloArray = find( collisionMatrix(:, startRadius + iRadiusTraveled) );
    
    % Loop through the non-empty velocity bins
    for jVelo = 1 : numel(jVeloArray)
       densityMedium = ( collisionMatrix(jVelo, ...
                                        startRadius + iRadiusTraveled) ...
                         / binVolumeArray(iRadiusTraveled) );
                     
       veloRatio = (iVelo - jVelo) / (iVelo + jVelo);
       
       colProb(iRadiusTraveled) = colProb(iRadiusTraveled) + ...
                                    (densityMedium * veloRatio);
    end
    
end

colProb = colProb .* (radiusDelta * scatterCS);

for i = 1 : numel(colProb)
    if colProb(i)
end

end

