function nParticleVeloDist = initialVelocityDistribution( ...
    uc, ucRatio, atomUC, nUCAblated, velo, initVeloDistWidth, ...
    absorbedRatio, energyLaser, heatTarget, energyExcitation, ...
    plotBool, nMin, nParticleTotal, distributionType )
%INITIALVELOCITYDISTRIBUTION Generate initial particle velocity ditribution
%   uc                  : Array of type UnitCells, containing unit cells defined
%                           in the class.
%   ucRatio             : Array of numbers same length as uc, containing the
%                           ratio between the unit cells in the target.
%   atomUC              : Array containing the unique elements present in the target
%   nUCAblated          : Array containing the number of unit cells ablated per
%                           unit cell present in the target.
%   velo                : Velocity array
%   initVeloDistWidth   : Velocity distribution width, only valid for
%                           normal and log-normal distributions.
%   absorbedRatio       : Ratio of laser energy absorbed by the target.
%   energyLaser         : Laser energy [J]
%   heatTarget          : Heat dissipation into the target [J]
%   energyExcitation    : Excitation energy per atom in the target [J]
%   plotBool            : Boolean that determined is results get plotted or not.
%   nMin                : Number of plasma particle threshold.
%   nParticleTotal      : Total number of ablated particles.
%   distributionType    : Initial velocity distribution type:
%                             'normal', 'log-normal', 'Maxwell-Boltzmann'

% Throw error if invalid distribution type is entered
if ~strcmp(distributionType, 'normal') && ~strcmp(distributionType, 'log-normal') ...
        && ~strcmp(distributionType, 'Maxwell-Boltzmann')
    error(['Invalid velocity distribution, choose: ' ...
           'normal, log-normal, or Maxwell-Boltzmann']);
end

% Initialize particle velocity distribution
nParticleVeloDist = zeros(numel(atomUC), numel(velo));

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
        veloRMS = sqrt( 2 * energyKinetic / uc(iUC).ELEMENTS(iElement).MASS );
        
        % Normal distribution
        if strcmp(distributionType, 'normal')
            veloDistTemp = exp(- (velo - veloRMS).^2 ./ (2 * initVeloDistWidth^2) ) ...
                                ./ ( sqrt(2*pi) * initVeloDistWidth );
        end
        
        % Log-normal distribution
        if strcmp(distributionType, 'log-normal')
            % Log-normal distribution mean
            mu = log( veloRMS^2 / sqrt(veloRMS^2 + initVeloDistWidth^2) );

            % Log-normal distribution variance
            sigma = sqrt( log( 1 + (initVeloDistWidth^2 / veloRMS^2) ) );

            veloDistTemp = exp( -( log(velo(2:end)) - mu ).^2 ./ (2 * sigma^2) ) ...
                ./ ( velo(2:end) .* (sqrt(2*pi) * sigma) );
        end

        % Maxwell-Boltzmann distribution
        if strcmp(distributionType, 'Maxwell-Boltzmann')
            veloDistTemp = velo.^2 .* exp(- ( (3*uc(iUC).ELEMENTS(iElement).MASS).*(velo.^2) ) ...
                                   ./ (4 * energyKinetic) );
        end
                        
        % Normalize by the number of ablated particles
        veloDistTemp = ( veloDistTemp ./ sum(veloDistTemp) ) ...
                       .* ( nUCAblated(iUC) * uc(iUC).AMOUNT(iElement) );
        
        % Add distribution to element
        nParticleVeloDist(elementIndex, :) = ...
            nParticleVeloDist(elementIndex, :) + veloDistTemp;
    end
end

% Normalize by total number of ablated particles
nParticleVeloDist = ( nParticleVeloDist ./ sum(nParticleVeloDist, 'all') ) ...
                    .* nParticleTotal;

% If plotting is enabled
if plotBool
    
    % Maximum required velocity
    iVeloMax = find( any(nParticleVeloDist > nMin, 1) , 1, 'last');
    
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
    fig = figure(numel(atomUC)+2);
    fig.Name = 'Initial velocity distribution';
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
    xlim([0 velo(iVeloMax)]);
end

end

