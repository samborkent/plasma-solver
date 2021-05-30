function nParticleVeloDist = initialVelocityDistribution( ...
    uc, ucRatio, atomUC, ablationVolume, velo, initVeloDistWidth, ...
    absorbedRatio, energyLaser, heatTarget, energyExcitation, plotInitVeloDistBool )
% function [nParticleVeloInit, velocityAverage] = initialVelocityDistribution( ...
%     plotBool, saveBool, velo, velDisWidth, nUCAblated, uc, ucRatio, ...
%     energyLaser, heatTarget, absorption, nParticleAngle, densityRatio )
%INITIALVELOCITYDISTRIBUTION Generate initial particle velocity ditribution

% Initialize particle velocity distribution
nParticleVeloDist = zeros(numel(atomUC), numel(velo));

% Number of ablated unit cells
nUCAblated = ucRatio .* ( ablationVolume / sum(ucRatio .* [uc.VOLUME]) );

% Loop through unit cells in the target
for iUC = 1 : numel(uc)   
    % Calculate the kinetic energy per atom in the unit cell
    energyKinetic = ( (absorbedRatio * energyLaser - heatTarget) ...
                      / ( nUCAblated(iUC) * sum(uc(iUC).AMOUNT) ) ) ...
                    - energyExcitation - uc(iUC).ENERGY_FORMATION;
    
    % Loop through the elements in the unit cell
    for iElement = 1 : numel(uc(iUC).AMOUNT)
        % Find index of element
        elementIndex = find(uc(iUC).ELEMENTS(iElement).NUMBER == [atomUC.NUMBER]);
        
        % Calculate the average velocity of each element in the unit cell
        veloAverage = sqrt( 2 * energyKinetic / uc(iUC).ELEMENTS(iElement).MASS);
        
        % Log-normal distribution mean
        mu = log( veloAverage^2 / sqrt(veloAverage^2 + initVeloDistWidth^2) );
        
        % Log-normal distribution variance
        sigma = sqrt( log( 1 + (initVeloDistWidth^2 / veloAverage^2) ) );
        
        % Initial particle velocity log- normal distribution
%         veloDistTemp = exp(- (velo - veloAverage).^2 ./ (2 * initVeloDistWidth^2) ) ...
%                             ./ ( sqrt(2*pi) * initVeloDistWidth );
        veloDistTemp = exp( -( log(velo(2:end)) - mu ).^2 ./ (2 * sigma^2) ) ...
            ./ ( velo(2:end) .* (sqrt(2*pi) * sigma) );
                        
        % Normalize by the number of ablated particles
        veloDistTemp = ( veloDistTemp ./ sum(veloDistTemp) ) ...
                       .* ( nUCAblated(iUC) * uc(iUC).AMOUNT(iElement) );
        
        % Add distribution to element
        nParticleVeloDist(elementIndex, 2:end) = ...
            nParticleVeloDist(elementIndex, 2:end) + veloDistTemp;
    end
end

% If plotting is enabled
if plotInitVeloDistBool
    
    % String for material formula
    formulaString = [];
    for i = 1 : numel(uc)
        if i == numel(uc)
            formulaString = [formulaString uc(i).FORMULA];
        else
            formulaString = [formulaString uc(i).FORMULA ' - '];
        end
    end
    
    % String for ratio of unit cells
    if numel(ucRatio) > 1
        ratioString = [];
        ucRatioTemp = ucRatio ./ min(ucRatio);
        for i = 1 : numel(ucRatio)
            if i == numel(uc)
                ratioString = [ratioString num2str(ucRatioTemp(i))];
            else
                ratioString = [ratioString num2str(ucRatioTemp(i)) ' : '];
            end
        end
    end

    % Initialize figure
    figure(numel(atomUC)+2);
    hold on;
    xlabel('Velocity [m/s]');
    ylabel('Number of particles');
    if numel(ucRatio) > 1
        title([ 'Initial velocity distribution of ' ...
                formulaString ' ( ' ratioString ' )' ]);
    else
        title([ 'Initial velocity distribution of ' ...
                formulaString ]);
    end
    
    % Loop through unique elements
    for iElement = 1 : numel(atomUC)
        % Plot result
        plot(velo, nParticleVeloDist(iElement, :), ...
             'LineWidth', 3, ...
             'DisplayName', atomUC(iElement).SYMBOL );
    end
    
    % Enable legend
    legend;
end

% % Make plot
% if plotBool
%     figVelDisInit = figure('Name', 'Initial Velocity Distribution');
%     hold on;
%     
%     if densityRatio == 1
%         title(['Single crystal ' formulaString]);
%     else
%         if numel(ucRatio) > 1
%             title(['Mixture ' formulaString ' (' ratioString ') ' ...
%                 '(Density: ' num2str(100 * densityRatio, 3) '%)']);
%         else
%             title(['Mixture ' formulaString ...
%                 ' (Density: ' num2str(100 * densityRatio, 3) '%)']);
%         end
%     end
%     
%     xlabel('Velocity [m/s]');
%     ylabel('Number of particles');
% %     xlim([0 max(velocity)]);
% %     ylim([0 2.5E13]);
% end
% 
% % Calculate constants
% nVelo = numel(velo);
% 
% atoms = [];
% 
% for i = 1 : numel(uc)
%     
%     % Pre-allocate
%     nParticleVeloInit = zeros( nVelo, numel(uc(i).AMOUNT) );
% 
%     % Initial kinetic energy of atoms (excitation is neglected) (eq. 1)
% %     energyKinetic = ((((energyLaser * absorption) - heatTarget) ...
% %                     / nUCAblated) - energyBinding) / nAtomUCSum;
%     energyKinetic = ((((energyLaser * absorption) - heatTarget) ...
%                     / nUCAblated(i)) - uc(i).ENERGY_FORMATION) ...
%                     / sum(uc(i).AMOUNT);
% 
%     for atom = numel(uc(i).AMOUNT) : -1 : 1
%         if (uc(i).AMOUNT(atom) > 0)           
%             % Initial average velocity
% %             velocityAverage = sqrt(2 * energyKinetic / atomUC(atom).MASS);
%             velocityAverage = sqrt(2 * energyKinetic / uc(i).ELEMENTS(atom).MASS);
%     %         velocityAverage = sqrt( ( (energyLaser * 0.9 * adsorptionC) / 3 / ...
%     %             nUCAblated - energyBinding) / 2 / atomUC(atom).MASS);
% 
%             for iVelo = 1 : nVelo
%                 % Initial particle distibution per velocity (eq. 3)
%                 nParticleVeloInit(iVelo, atom) = exp(-(velo(iVelo) - velocityAverage)^2 / ...
%                     (2 * velDisWidth^2)) / (velDisWidth * sqrt(2 * pi));
%             end % For velo
% 
%             % Normalization factor
% %             nNorm = (nParticleAngle * (uc(i).AMOUNT(atom) / sum(uc(i).AMOUNT))) ...
% %                         / sum(nParticleVeloInit(:, atom));
%             nNorm = ( nParticleAngle * (ucRatio(i) / sum(ucRatio)) ...
%                     * (uc(i).AMOUNT(atom) / sum(uc(i).AMOUNT)) )...
%                     / sum(nParticleVeloInit(:, atom));
% 
%             % Normalize velocity distribution
%             nParticleVeloInit(:, atom) = nParticleVeloInit(:, atom) .* nNorm;
% 
%             % Plot results
%             if plotBool
%                 plot(velo, nParticleVeloInit(:, atom), ...
%                     'DisplayName', uc(i).ELEMENTS(atom).SYMBOL);
% %                 bar(velo, nParticleVeloInit(:, atom), 'stacked', ...
% %                     'LineStyle', 'none', ...
% %                     'DisplayName', uc(i).ELEMENTS(atom).SYMBOL);
%             end % If plotBool
%         end % If there are any atoms
%     end % For atoms in target component
%     
% end % For target component
% 
% % Generate legend
% if plotBool
%     legend;
%     set(findobj('Type','line'), 'LineWidth', 3)
%     
%     % Save figure
%     if saveBool
%         formulaString(formulaString == '_') = [];
%         ratioString(ratioString == ':') = '-';
%         
%         if numel(ucRatio) > 1
%             saveas(figVelDisInit, [pwd '/results/' formulaString '_' ...
%                 ratioString '_' num2str(densityRatio, 3) '.png']);
%         else
%             if densityRatio == 1
%                 saveas(figVelDisInit, [pwd '/results/' formulaString ...
%                     '.png']);
%             else
%                 saveas(figVelDisInit, [pwd '/results/' formulaString '_' ...
%                     '_' num2str(densityRatio, 3) '.png']);
%             end
%         end
%     end
% end

end

