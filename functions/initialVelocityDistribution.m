function [nParticleVeloInit, velocityAverage] = initialVelocityDistribution( ...
    plotBool, saveBool, velo, velDisWidth, nUCAblated, uc, ucRatio, ...
    energyLaser, heatTarget, absorption, nParticleAngle, densityRatio )
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

formulaString = [];

for i = 1 : numel(uc)
    if i == numel(uc)
        formulaString = [formulaString uc(i).FORMULA];
    else
        formulaString = [formulaString uc(i).FORMULA ' - '];
    end
end

ratioString = [];

for i = 1 : numel(ucRatio)
    if i == numel(uc)
        ratioString = [ratioString num2str(ucRatio(i))];
    else
        ratioString = [ratioString num2str(ucRatio(i)) ':'];
    end
end

% Make plot
if plotBool
    figVelDisInit = figure('Name', 'Initial Velocity Distribution');
    hold on;
    
    if densityRatio == 1
        title(['Single crystal ' formulaString]);
    else
        if numel(ucRatio) > 1
            title(['Mixture ' formulaString ' (' ratioString ') ' ...
                '(Density: ' num2str(100 * densityRatio, 3) '%)']);
        else
            title(['Mixture ' formulaString ...
                ' (Density: ' num2str(100 * densityRatio, 3) '%)']);
        end
    end
    
    xlabel('Velocity [m/s]');
    ylabel('Number of particles');
%     xlim([0 max(velocity)]);
%     ylim([0 2.5E13]);
end

% Calculate constants
nVelo = numel(velo);

atoms = [];

for i = 1 : numel(uc)
    
    % Pre-allocate
    nParticleVeloInit = zeros( nVelo, numel(uc(i).AMOUNT) );

    % Initial kinetic energy of atoms (excitation is neglected) (eq. 1)
%     energyKinetic = ((((energyLaser * absorption) - heatTarget) ...
%                     / nUCAblated) - energyBinding) / nAtomUCSum;
    energyKinetic = ((((energyLaser * absorption) - heatTarget) ...
                    / nUCAblated(i)) - uc(i).ENERGY_FORMATION) ...
                    / sum(uc(i).AMOUNT);

    for atom = numel(uc(i).AMOUNT) : -1 : 1
        if (uc(i).AMOUNT(atom) > 0)           
            % Initial average velocity
%             velocityAverage = sqrt(2 * energyKinetic / atomUC(atom).MASS);
            velocityAverage = sqrt(2 * energyKinetic / uc(i).ELEMENTS(atom).MASS);
    %         velocityAverage = sqrt( ( (energyLaser * 0.9 * adsorptionC) / 3 / ...
    %             nUCAblated - energyBinding) / 2 / atomUC(atom).MASS);

            for iVelo = 1 : nVelo
                % Initial particle distibution per velocity (eq. 3)
                nParticleVeloInit(iVelo, atom) = exp(-(velo(iVelo) - velocityAverage)^2 / ...
                    (2 * velDisWidth^2)) / (velDisWidth * sqrt(2 * pi));
            end % For velo

            % Normalization factor
%             nNorm = (nParticleAngle * (uc(i).AMOUNT(atom) / sum(uc(i).AMOUNT))) ...
%                         / sum(nParticleVeloInit(:, atom));
            nNorm = ( nParticleAngle * (ucRatio(i) / sum(ucRatio)) ...
                    * (uc(i).AMOUNT(atom) / sum(uc(i).AMOUNT)) )...
                    / sum(nParticleVeloInit(:, atom));

            % Normalize velocity distribution
            nParticleVeloInit(:, atom) = nParticleVeloInit(:, atom) .* nNorm;

            % Plot results
            if plotBool
                plot(velo, nParticleVeloInit(:, atom), ...
                    'DisplayName', uc(i).ELEMENTS(atom).SYMBOL);
%                 bar(velo, nParticleVeloInit(:, atom), 'stacked', ...
%                     'LineStyle', 'none', ...
%                     'DisplayName', uc(i).ELEMENTS(atom).SYMBOL);
            end % If plotBool
        end % If there are any atoms
    end % For atoms in target component
    
end % For target component

% Generate legend
if plotBool
    legend;
    set(findobj('Type','line'), 'LineWidth', 3)
    
    % Save figure
    if saveBool
        formulaString(formulaString == '_') = [];
        ratioString(ratioString == ':') = '-';
        
        if numel(ucRatio) > 1
            saveas(figVelDisInit, [pwd '/results/' formulaString '_' ...
                ratioString '_' num2str(densityRatio, 3) '.png']);
        else
            if densityRatio == 1
                saveas(figVelDisInit, [pwd '/results/' formulaString ...
                    '.png']);
            else
                saveas(figVelDisInit, [pwd '/results/' formulaString '_' ...
                    '_' num2str(densityRatio, 3) '.png']);
            end
        end
    end
end

end

