function [nVelDisInit, figVelDisInit] = initialVelocityDistribution( ...
    vel, velDisWidth, nUCAblated, nAtomsUC, massSpecies, energyLaser, ...
    energyFormation, energyExcitation, heatTarget, adsorptionC )
%INITIALVELOCITYDISTRIBUTION Generate initial particle velocity ditribution
%   vel                 Velocity array
%                           Example: [0 : (VEL_MIN / VEL_DELTA) : V_MAX]
%   velDisWidth         Width of velocity ditribution
%   nUCAblated          Number of ablated unit cells per laser pulse
%   nAtomUC             Array with amount of each atom in the unit cell
%                           Example: [2 1 3] for Li2TiO3
%   massAtoms           Array with the mass of the atoms in the unit cell
%   energyLaser         Laser energy [J]
%   energyFormation     Formation energy of unit cell [J]
%   energyExcitation    Excitation energy array of atoms in unit cell [J]
%   heatTarget          Heat loss into the target [J]
%   adsorptionC         Adsorption coefficient of the target [0 : 1]
%


figVelDisInit = figure;
hold on;

% Calculate constants
N_VEL = numel(vel);

% Pre-allocate
nVelDisInit = zeros(1, N_VEL);

for atom = 1 : numel(nAtomsUC)
    % Initial kinetic energy of species (eq. 1)
    energyKinetic = (((energyLaser / adsorptionC) - ...
        heatTarget) / nUCAblated) - energyFormation - ...
        energyExcitation(atom);
    
    % Initial average velocity
    velAvgInit = sqrt(2 * energyKinetic / massSpecies(atom));
    
    for i = 1 : N_VEL
        % Initial particle distibution per velocity (eq. 3)
        nVelDisInit(i) = exp(-(vel(i) - velAvgInit)^2 / ...
            (2 * velDisWidth^2)) / (velDisWidth * sqrt(2 * pi));
    end
    
    % Normalization factor
    nNorm = (nUCAblated * nAtomsUC(atom)) / sum(nVelDisInit);
    
    % Normalize velocity distribution
    nVelDisInit = nVelDisInit .* nNorm;
    
    plot(vel, nVelDisInit);
end

end

