classdef PhysicalConstants
    %PHYSICALCONSTANTS Class containing physical constants
    
    properties (Constant)
        %% Nature constants
        BOLTZMANN = 1.38064852*10^-23;  % Boltzman constant [m2 kg / s2 K]
        PLANCK = 6.62607004E-34;        % Planck constant [m2 kg / s]
        SOL = 299792458;                % Speed of light [m / s]
        
        %% SI conversion constants
        AMU = 1.66053906660E-27;        % Atomic mass unit [kg]
        EV  = 1.602176634E-19;          % Electron volt [J]
        HA  = 4.3597447222071E-18;      % Hatree [J]
      
        %% Unit cell volumes [m]
        % TiO2, tetragonal => a^2 * c
        % Source: "Numerical modelling of the plasma plume propagation ...
        %   and oxidation during pulsed laser deposition of complex ...
        %   oxide thin films.", 2020, T. Wijnands, et al.
%         UC_VOL_TiO2      = (3.78E-10)^2 * 9.51E-10;   % Paper
        UC_VOL_TiO2 = (4.5937E-10^2) * 2.9587E-10;      % Script
        
        % SrTiO3, cubic => a^3
        UC_VOL_SrTiO3    = (3.90E-10)^3;
        
        % Li4Ti5O12, cubic => a^3
        % Source: "Doubling Reversible Capacities in Epitaxial Li4Ti5O12
        %   Thin Film Anodes for Microbatteries", 2019, D.M. Cunha et al.
        UC_VOL_Li4Ti5O12 = (8.36E-10)^3;
        
        %% Unit cell formation energies [J]
        % Source: https://materialsproject.org/materials/mp-....../
        
        % Spinel Li4Ti5O12
        % Source: 'Doubling Reversible Capacities in Epitaxial ...
        %   Li4Ti5O12 Thin Film Anodes for Microbatteries', 2019, ...
        %   D.M. Cunha et al.
        UC_EF_Li4Ti5O12 = 3.2 * PhysicalConstants.EV;
    end
end
