classdef UnitCells
    %UNITCELLS Class holding material unit cell data
    % Source: materialsproject.org
    
    properties
        TiO2 = Material( ...                                         % mp-2657
            'TiO_2', ...                                             % Formula
            [PeriodicTable.Ti PeriodicTable.O], ...                  % Elements
            [1 2], ...                                               % Amount of each atom
            3.475 * PhysicalConstants.EV, ...                        % Formation energy [J]
            4.653E-10, ...                                           % Lattice constant a [m]
            4.653E-10, ...                                           % Lattice constant b [m]
            2.969E-10, ...                                           % Lattice constant c [m]
            4130 ...                                                 % Single crystal density [kg / m^3]
           );
        SrTiO3 = Material( ...                                       % mp-5229
            'SrTiO_3', ...                                           % Formula
            [PeriodicTable.Sr PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [1 1 3], ...                                             % Amount of each atom
            3.561 * PhysicalConstants.EV, ...                        % Formation energy [J]
            3.945E-10, ...                                           % Lattice constant a [m]
            3.945E-10, ...                                           % Lattice constant b [m]
            3.945E-10, ...                                           % Lattice constant c [m]
            4960 ...                                                 % Single crystal density [kg / m^3]
           );
        Li2O = Material( ...                                         % mp-1960
            'Li_2O', ...                                             % Formula
            [PeriodicTable.Li PeriodicTable.O], ...                  % Elements
            [2 1], ...                                               % Amount of each atom
            2.067 * PhysicalConstants.EV, ...                        % Formation energy [J]
            3.294E-10, ...                                           % Lattice constant a [m]
            3.294E-10, ...                                           % Lattice constant b [m]
            3.294E-10, ...                                           % Lattice constant c [m]
            1960 ...                                                 % Single crystal density [kg / m^3]
           );   
        LiO2 = Material( ...                                         % mp-1018789
            'LiO_2', ...                                             % Formula
            [PeriodicTable.Li PeriodicTable.O], ...                  % Elements
            [1 2], ...                                               % Amount of each atom
            1.033 * PhysicalConstants.EV, ...                        % Formation energy [J]
            2.961E-10, ...                                           % Lattice constant a [m]
            3.986E-10, ...                                           % Lattice constant b [m]
            4.896E-10, ...                                           % Lattice constant c [m]
            2240 ...                                                 % Single crystal density [kg / m^3]
           );       
        LiTi2O4 = Material( ...                                      % mp-5670
            'LiTi_2O_4', ...                                         % Formula
            [PeriodicTable.Li PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [1 2 4], ...                                             % Amount of each atom
            3.255 * PhysicalConstants.EV, ...                        % Formation energy [J]
            5.986E-10, ...                                           % Lattice constant a [m]
            5.986E-10, ...                                           % Lattice constant b [m]
            5.986E-10, ...                                           % Lattice constant c [m]
            3650 ...                                                 % Single crystal density [kg / m^3]
           );
        Li2TiO3 = Material( ...
            'Li_2TiO_3', ...                                         % Formula
            [PeriodicTable.Li PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [2 1 3], ...                                             % Amount of each atom
            2.979 * PhysicalConstants.EV, ...                        % Formation energy [J]
            5.0597E-10, ...                                          % Lattice constant a [m] hindawi.com/journals/jchem/2019/1727859/
            5.0597E-10, ...                                          % Lattice constant b [m]
            5.087E-10, ...                                           % Lattice constant c [m]
            3420 ...                                                 % Single crystal density [kg / m^3]
           );
        Li4Ti5O12 = Material( ...
            'Li_4Ti_5O_1_2', ...                                     % Formula
            [PeriodicTable.Li PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [4 5 12], ...                                            % Amount of each atom
            3.192 * PhysicalConstants.EV, ...                        % Formation energy [J]
            8.3557E-10, ...                                          % Lattice constant a [m] hindawi.com/journals/jchem/2019/1727859/
            8.3557E-10, ...                                          % Lattice constant b [m]
            8.3557E-10, ...                                          % Lattice constant c [m]
            3290 ...                                                 % Single crystal density [kg / m^3]
           );
    end
end

