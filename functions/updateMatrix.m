function particleMatrix = updateMatrix( particleMatrix, nVelo, nSpecies, kMax, ...
                                        preserveParticlesBool )
%UPDATEMATRIX Update the non-colliding particle positions.
%   If preserveParticlesBool is true, particles will remain in the last
%   radial bin and the total number of particles in the particle matrix
%   will be conserved. Has a significant hit on performance.
%   If preserveParticlesBool is false, particles will disappear after
%   reaching the last radial bin.

for iVelo = nVelo : -1 : 2
%% Calculations per plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles
% Skip zero velocity

    if preserveParticlesBool
        % Sum all particles that will pass the last radial bin and insert
        %   change their position so they will remain in the last radial bin.
        particleMatrix(:, :, iVelo, end-iVelo+1) = ...
            sum( particleMatrix(:, :, iVelo, end-iVelo+1:end), 'all' );
    end
    
    % Move particles to new radial bin based on position
    particleMatrix(:, :, iVelo, iVelo:end) = ...
        particleMatrix(:, :, iVelo, 1:end-iVelo+1);
    
    % Remove particles from previous position
    particleMatrix(:, :, iVelo, 1:iVelo-1) = zeros(nSpecies, kMax, 1, iVelo-1);
end
    
end

