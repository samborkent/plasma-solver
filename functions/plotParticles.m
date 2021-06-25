function plotParticles( particleMatrix, plotTimes, iTime, time, atomUC, ...
    nMin, radius, binVolume, smoothPlotBool, plotDensityBool, plotLogBool )
%PLOTRESULT Plot the particle propagation.

% If the current time is one of the times to plot
if any(iTime == plotTimes)
    
% Indices for number of collisions
%   (1) : Non-collided
%   (2) : Collided
%   (3) : Sum

% Array holding line colors for each number of collisions
colorArray = ['b' 'r' 'k'];

% Array holding display name for each number of collisions
plotArray = {'Uncollided', 'Collided', 'Total'};

% Get dimensional sizes to prevent having to pass many global variables
%   as arguments
nSpecies = size(particleMatrix, 1);
nK = size(particleMatrix, 2);
nRadius = size(particleMatrix, 4);

% Get subplot index
plotTimeIndex = find( iTime == plotTimes );

% Loop through plasma species
for iSpecies = 1 : nSpecies
    % Initialize figure
    fig = figure(iSpecies);
    fig.Position = [250 100 1247 288];

    % Set figure name
    if iSpecies == 1
        fig.Name = 'Background propagation';
    else
        fig.Name = [atomUC(iSpecies-1).SYMBOL ' propagation'];
    end

    % Select subplot
    subplot(1, numel(plotTimes), plotTimeIndex);
    hold on;

    % Set title
    title([num2str(time(iTime), 3) ' s']);

    % Initialize sum of all collisions array
    nParticleRadiusSum = zeros(1, nRadius);

    % Loop through number of collisions
    for k = 1 : nK
        % Calculate total number of particles per radial bin
        nParticleRadius = sum( squeeze( particleMatrix(iSpecies, k, :, :) ), 1 );

        % If plot smoothing is enabled and the density of the
        %   background gas is not plotted
        if smoothPlotBool && ~(plotDensityBool && (iSpecies == 1))
            % Smooth data
            nParticleRadius = ...
                smoothdata( nParticleRadius, ...            % Data
                            2, ...                          % Dimension
                            'gaussian', ...                 % Type
                            round(nRadius * max(radius)) ); % Smoothing width
        end

        % Add to sum of all collisions array
        nParticleRadiusSum = nParticleRadiusSum + nParticleRadius;

        % If background gas species and density plot is enabled
        if (iSpecies == 1) && plotDensityBool
            plot( radius(nParticleRadius >= nMin), nParticleRadius(nParticleRadius >= nMin) ./ binVolume(nParticleRadius >= nMin), ...
                  'Color', colorArray(k), ...
                  'LineWidth', 2, ...
                  'DisplayName', plotArray{k} );
        else
            plot( radius(nParticleRadius >= nMin), nParticleRadius(nParticleRadius >= nMin), ...
                  'Color', colorArray(k), ...
                  'LineWidth', 2, ...
                  'DisplayName', plotArray{k} );
        end
    end % Number of collisions loop

    % If background gas species and density plot is enabled
    if (iSpecies == 1) && plotDensityBool
        plot( radius(nParticleRadiusSum >= nMin), nParticleRadiusSum(nParticleRadiusSum >= nMin) ./ binVolume(nParticleRadiusSum >= nMin), ...
              'Color', colorArray(end), ...
              'LineWidth', 2, ...
              'DisplayName', plotArray{end} );
    else
        plot( radius(nParticleRadiusSum >= nMin), nParticleRadiusSum(nParticleRadiusSum >= nMin), ...
              'Color', colorArray(end), ...
              'LineWidth', 2, ...
              'DisplayName', plotArray{end} );
    end

    % Set log-scale if selected
    if plotLogBool
        set(gca, 'YScale', 'log');
    end

    % Enable legend for last subplot
    if plotTimeIndex == numel(plotTimes)
        legend;
        legend('boxoff');
        legend('Location', 'northeast');
    end
    
    % Set limits of axis equal for all subplots
    linkaxes(findall(figure(iSpecies), 'type', 'axes'))
    
    % Set x-limits
    xlim([0 max(radius)]);

    hold off;

end % Species loop

end % If plotTimes

end

