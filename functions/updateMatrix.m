function particleMatrix = updateMatrix( particleMatrix, ...
                                        nVelo, iVeloZero, ...
                                        keepParticlesBool )
%UPDATEMATRIX Update the non-colliding particle positions.
%   If preserveParticlesBool is true, particles will remain in the last
%   radial bin and the total number of particles in the particle matrix
%   will be conserved. Has a significant hit on performance.
%   If preserveParticlesBool is false, particles will disappear after
%   exceeding the first or last radial bin.

for iVelo = nVelo : -1 : 1
%% Calculations per plasma velocity bin
% Loop backwards as fast particles travel in front of slower particles

    % Skip zero velocity
    if iVelo == iVeloZero
        continue
    end
    
    % If number of particles is conserved
    if keepParticlesBool
        % Positive velocity
        if iVelo > iVeloZero
            % Add all disappearing particles to the last radial bin with
            %   zero velocity
            particleMatrix(:, :, iVeloZero, end) = ...
                particleMatrix(:, :, iVeloZero, end) ...
                + sum( particleMatrix(:, :, iVelo, end-(iVelo-iVeloZero)+1:end), 4 );
        % Negative velocity
        elseif iVelo < iVeloZero
            % Add all disappearing partilces to the first radial bin with
            %   zero veloctiy
            particleMatrix(:, :, iVeloZero, 1) = ...
                particleMatrix(:, :, iVeloZero, 1) ...
                + sum( particleMatrix(:, :, iVelo, 1:iVeloZero-iVelo), 4 );
        end
    end

    % Move particles to new position based on velocity
    particleMatrix(:, :, iVelo, :) = ...
        circshift( particleMatrix(:, :, iVelo, :), iVelo - iVeloZero, 4 );
    
    % Remove particles from previous position
    if iVelo > iVeloZero
        particleMatrix(:, :, iVelo, 1:iVelo-iVeloZero) = 0;
    elseif iVelo < iVeloZero
        particleMatrix(:, :, iVelo, end-(iVeloZero-iVelo)+1:end) = 0;
    end

end
    
end

