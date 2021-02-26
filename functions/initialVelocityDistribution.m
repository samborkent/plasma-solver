function nParticleVeloInit = initialVelocityDistribution( ...
    plotBool, saveBool, velo, velDisWidth, nUCAblated, uc, ...
    energyLaser, energyBinding, heatTarget, absorption, nParticleAngle )
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
%   energyBinding       Binding energy of unit cell [J]
%   heatTarget          Heat loss into the target [J]
%   absorption          Absorption coefficient of the target [0 : 1]
%

% Make plot
if plotBool
    figVelDisInit = figure;
    hold on;
    title(uc.FORMULA);
    xlabel('Velocity [m/s]');
    ylabel('Number of particles');
%     xlim([0 max(velocity)]);
%     ylim([0 2.5E13]);
end

% Calculate constants
nVelo = numel(velo);

% Pre-allocate
nParticleVeloInit = zeros(nVelo, numel(uc.ELEMENTS));

% Initial kinetic energy of atoms (excitation is neglected) (eq. 1)
energyKinetic = ((((energyLaser * absorption) - heatTarget) ...
                / nUCAblated) - energyBinding) / sum(uc.AMOUNT);

for atom = numel(uc.ELEMENTS) : -1 : 1
    if (uc.AMOUNT(atom) > 0)       
        % Initial average velocity
        velocityAverage = sqrt(2 * energyKinetic / uc.ELEMENTS(atom).MASS);
%         velocityAverage = sqrt( ( (energyLaser * 0.9 * adsorptionC) / 3 / ...
%             nUCAblated - energyBinding) / 2 / atomUC(atom).MASS);

        for iVelo = 1 : nVelo
            % Initial particle distibution per velocity (eq. 3)
            nParticleVeloInit(iVelo, atom) = exp(-(velo(iVelo) - velocityAverage)^2 / ...
                (2 * velDisWidth^2)) / (velDisWidth * sqrt(2 * pi));
        end

        % Normalization factor
        nNorm = (nParticleAngle * (uc.AMOUNT(atom) / sum(uc.AMOUNT))) ...
                    / sum(nParticleVeloInit(:, atom));

        % Normalize velocity distribution
        nParticleVeloInit(:, atom) = nParticleVeloInit(:, atom) .* nNorm;

        % Plot results
        if plotBool
%             plot(velocity, velocityDistributionInitial, 'DisplayName', atomUC(atom).SYMBOL);
            bar(velo, nParticleVeloInit(:, atom), 'grouped', ...
                'LineStyle', 'none', 'DisplayName', uc.ELEMENTS(atom).SYMBOL);
        end
    end
end

% Generate legend
if plotBool
    legend;
    set(findobj('Type','line'), 'LineWidth', 3)
    
    % Save figure
    if saveBool
        saveas(figVelDisInit, [pwd '/results/' uc.FORMULA '.png']);
    end
end

end

