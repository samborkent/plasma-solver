classdef UnitCells
    %UNITCELLS Class holding material unit cell data
    % Source: materialsproject.org
    
    properties
        TiO2 = Material( ...                                         % mp-2657
            'TiO_2', ...                                             % Formula
            [PeriodicTable.Ti PeriodicTable.O], ...                  % Elements
            [2 4], ...                                               % Amount of each atom
            66.292E-30, ...                                          % Unit cell volume [m^3]
            4130, ...                                                % Single crystal density [kg / m^3]
            3.464 * PhysicalConstants.EV, ...                        % Formation energy per atom [J]
            56.3693 * PhysicalConstants.EV ...                       % Total energy unit cell [J]
           );
        SrTiO3 = Material( ...                                       % mp-5229
            'SrTiO_3', ...                                           % Formula
            [PeriodicTable.Sr PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [1 1 3], ...                                             % Amount of each atom
            61.402E-30, ...                                          % Unit cell volume [m^3]
            4960, ...                                                % Single crystal density [kg / m^3]
            3.551 * PhysicalConstants.EV, ...                        % Formation energy per atom [J]
            42.1857 * PhysicalConstants.EV ...                       % Total energy unit cell [J]
           );
        Li2O = Material( ...                                         % mp-1960
            'Li_2O', ...                                             % Formula
            [PeriodicTable.Li PeriodicTable.O], ...                  % Elements
            [2 1], ...                                               % Amount of each atom
            25.280E-30, ...                                          % Unit cell volume [m^3]
            1960, ...                                                % Single crystal density [kg / m^3]
            2.067 * PhysicalConstants.EV, ...                        % Formation energy per atom [J]
            14.9506 * PhysicalConstants.EV ...                       % Total energy unit cell [J]
           );   
        LiTi2O4 = Material( ...                                      % mp-5670
            'LiTi_2O_4', ...                                         % Formula
            [PeriodicTable.Li PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [2 4 8], ...                                             % Amount of each atom
            25.280E-30, ...                                          % Unit cell volume [m^3]
            3650, ...                                                % Single crystal density [kg / m^3]
            3.246 * PhysicalConstants.EV, ...                        % Formation energy per atom [J]
            120.4304 * PhysicalConstants.EV ...                      % Total energy unit cell [J]
           );
        Li4Ti5O12 = Material( ...                                    % mp-685194
            'Li_4Ti_5O_1_2', ...                                     % Formula
            [PeriodicTable.Li PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [8 10 24], ...                                           % Amount of each atom
            595.327E-30, ...                                         % Unit cell volume [m^3]
            3380, ...                                                % Single crystal density [kg / m^3]
            3.210 * PhysicalConstants.EV, ...                        % Formation energy per atom [J]
            347.7836 * PhysicalConstants.EV ...                      % Total energy unit cell [J]
           );
       LiLaTi2O6 = Material( ...                                     % mp-1222315
            'LiLaTi_2O_6', ...                                       % Formula
            [PeriodicTable.Li PeriodicTable.La PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [1 1 2 6], ...                                           % Amount of each atom
            117.729E-30, ...                                         % Unit cell volume [m^3]
            4760, ...                                                % Single crystal density [kg / m^3]
            3.383 * PhysicalConstants.EV, ...                        % Formation energy per atom [J]
            86.1521 * PhysicalConstants.EV ...                       % Total energy unit cell [J]
           );
       LiMn2O4 = Material( ...                                       % mp-22584
            'LiMn_2O_4', ...                                         % Formula
            [PeriodicTable.Li PeriodicTable.Mn PeriodicTable.O], ... % Elements
            [2 4 8], ...                                             % Amount of each atom
            149.369E-30, ...                                         % Unit cell volume [m^3]
            4020, ...                                                % Single crystal density [kg / m^3]
            2.023 * PhysicalConstants.EV, ...                        % Formation energy per atom [J]
            108.3734 * PhysicalConstants.EV ...                      % Total energy unit cell [J]
           );
    end
end

