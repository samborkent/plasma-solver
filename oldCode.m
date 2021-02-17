%% Clear workspace
clear
clearvars -except files nn
close all
clc

% Select execution time profiler
profile off
warning('OFF')

%--------------------------------------------------------------------------
% Plot settings
%--------------------------------------------------------------------------

ThreeD      = false;            % Plot 3D results: true / false
Times_3D    = [5, 15, 20, 30];  % Times for which to generate 3D plots [s]
Arr_sum    = false;            % Calculate arriving particles:  true / false
Plot_1D     = false;            % Plot 1D results of middle of the plume: true / false
Plot_BG     = false;            % Show bg in plot: true / false

%--------------------------------------------------------------------------
% Simulation setting
%--------------------------------------------------------------------------

% Restrictions (Check values)
enableOxidation = false;    %?% Enable particle oxidation (true = no oxidization)
enableDioxides  = false;    % Let particles form dioxides: true / false
nMinPerBin = 1;             % Minimal number of particles per bin
nMinPerAngle = 1E7;         % Minimal N per angle (angles with lower N are not computed)
nCollisionsMax = 7;         % Maximum number of collisions per particle

% Smoothing parameters
Smooth_Fraction = 1E4;  % Fraction to smooth % ???
Smooth_Value    = 15;   % Moving average smooth value (default: 15)