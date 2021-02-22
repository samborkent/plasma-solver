classdef UnitCells
    %UNITCELLS Class holding material unit cell data
    % Source: materialsproject.org
    
    properties
        TiO2 = Material( ...                                         % mp-2657
            'TiO2', ...                                              % Formula
            [PeriodicTable.Ti PeriodicTable.O], ...                  % Elements
            [1 2], ...                                               % Amount of each atom
            3.475 * PhysicalConstants.EV, ...                        % Formation energy [J]
            1.772 * PhysicalConstants.EV, ...                        % Band gap energy [J]
            4.653E-10, ...                                           % Lattice constant a [m]
            4.653E-10, ...                                           % Lattice constant b [m]
            2.969E-10, ...                                           % Lattice constant c [m]
            4130 ...                                                 % Density [kg / m^3]
           );
        SrTiO3 = Material( ...                                       % mp-5229
            'SrTiO3', ...                                            % Formula
            [PeriodicTable.Sr PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [1 1 3], ...                                             % Amount of each atom
            3.561 * PhysicalConstants.EV, ...                        % Formation energy [J]
            1.872 * PhysicalConstants.EV, ...                        % Band gap energy [J]
            3.945E-10, ...                                           % Lattice constant a [m]
            3.945E-10, ...                                           % Lattice constant b [m]
            3.945E-10, ...                                           % Lattice constant c [m]
            4960 ...                                                 % Density [kg / m^3]
           );
       Li2TiO3 = Material( ...
            'Li2TiO3', ...                                           % Formula
            [PeriodicTable.Li PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [2 1 3], ...                                             % Amount of each atom
            2.979 * PhysicalConstants.EV, ...                        % Formation energy [J]
            2.418 * PhysicalConstants.EV, ...                        % Band gap energy [J]
            5.0597E-10, ...                                          % Lattice constant a [m] hindawi.com/journals/jchem/2019/1727859/
            5.0597E-10, ...                                          % Lattice constant b [m]
            5.087E-10, ...                                           % Lattice constant c [m]
            3420 ...                                                 % Density [kg / m^3]
           );
       Li4Ti5O12 = Material( ...
            'Li4Ti5O12', ...                                         % Formula
            [PeriodicTable.Li PeriodicTable.Ti PeriodicTable.O], ... % Elements
            [4 5 12], ...                                            % Amount of each atom
            3.192 * PhysicalConstants.EV, ...                        % Formation energy [J]
            3.015 * PhysicalConstants.EV, ...                        % Band gap energy [J]
            8.3557E-10, ...                                          % Lattice constant a [m] hindawi.com/journals/jchem/2019/1727859/
            8.3557E-10, ...                                          % Lattice constant b [m]
            8.3557E-10, ...                                          % Lattice constant c [m]
            3290 ...                                                 % Density [kg / m^3]
           );
    end
end

