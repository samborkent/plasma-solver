function [plasmaTemp, bgTemp] = averageTemperature( particleMatrix, velo, mass )
%AVERAGETEMPERATURE Calculate the average temeparture of the plasma and the
%   background gas.

% Number of background particles per velocity
nBGVelo = squeeze(sum(particleMatrix(1, :, :, :), [2 4]))';

% Mean velocity of background gas
bgVelo = sum( nBGVelo .* velo ) / sum( nBGVelo );

% Mean kinetic energy per background gas particle
bgKinetic = 0.5 * mass(1) * bgVelo^2;

% Thermodynamic temperature of background gas
bgTemp = (2/3) * bgKinetic / PhysicalConstants.BOLTZMANN;

% Number of plasma particles per species
nPlasmaSpecies = zeros(1, size(particleMatrix, 1) - 1);

% Mean kinetic energy per plasma species
plasmaKinetic = zeros(1, size(particleMatrix, 1) - 1);

% Loop through plasma species
for iSpecies = 2 : size(particleMatrix, 1)
    % Number of plasma particles per velocity of this species
    nPlasmaVelo = squeeze(sum(particleMatrix(iSpecies, :, :, :), [2 4]))';
    
    % Store total number of particles of this species
    nPlasmaSpecies(iSpecies-1) = sum(nPlasmaVelo, 2);
    
    % Mean velocity of this species
    plasmaVelo = sum( nPlasmaVelo .* velo ) / nPlasmaSpecies(iSpecies-1);
    
    % Mean kinetic energy of this species
    plasmaKinetic(iSpecies-1) = 0.5 * mass(iSpecies) * plasmaVelo^2;
end

% Mean kinetic energy of all plasma particles
meanPlasmaKinetic = sum( plasmaKinetic .* nPlasmaSpecies ) / sum( nPlasmaSpecies );

% Thermodynamic temperature of plasma
plasmaTemp = (2/3) * meanPlasmaKinetic / PhysicalConstants.BOLTZMANN;

end