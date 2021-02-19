function [nVelDisInit, energyKinetic, figVelDisInit] = initialVelocityDistribution( ...
    plotBool, velocity, velDisWidth, nUCAblated, atomUC, nAtomUC, ...
    energyLaser, energyBinding, heatTarget, adsorptionC, nParticleAngle )
%INITIALVELOCITYDISTRIBUTION Generate initial particle velocity ditribution
%   plotBool            Boolean which determines if results will be plotted
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

% Make plot
if plotBool
    figVelDisInit = figure;
    hold on;
    xlabel('Velocity [m/s]');
    ylabel('Number of particles');
end

% Calculate constants
nVelocity = numel(velocity);

% Pre-allocate
nVelDisInit = zeros(1, nVelocity);

% Initial kinetic energy of atoms (eq. 1)
energyKinetic = (((energyLaser / adsorptionC) - heatTarget) ...
                / nUCAblated) - energyBinding - atomUC(1).ENERGY_FI ...
                - atomUC(2).ENERGY_FI - atomUC(3).ENERGY_FI;

for atom = 1 : numel(nAtomUC)
    if (nAtomUC(atom) > 0)       
        % Initial average velocity
        velocityAverage = sqrt(2 * energyKinetic / atomUC(atom).MASS);
%         velocityAverage = sqrt( ( (energyLaser * 0.9 * adsorptionC) / 3 / ...
%             nUCAblated - energyBinding) / 2 / atomUC(atom).MASS);

        for i = 1 : nVelocity
            % Initial particle distibution per velocity (eq. 3)
            nVelDisInit(i) = exp(-(velocity(i) - velocityAverage)^2 / ...
                (2 * velDisWidth^2)) / (velDisWidth * sqrt(2 * pi));
        end

        % Normalization factor
    %     nNorm = (nParticleAngle * nAtomUC(atom)) / sum(nVelDisInit);
        nNorm = nParticleAngle * (nAtomUC(atom) / sum(nAtomUC)) / sum(nVelDisInit);

        % Normalize velocity distribution
        nVelDisInit = nVelDisInit .* nNorm;

        % Plot results
        if plotBool
            plot(velocity, nVelDisInit, 'DisplayName', atomUC(atom).SYMBOL);
        end
    end
end

% Generate legend
if plotBool
    legend;
    set(findobj('Type','line'), 'LineWidth', 3)
end

end

