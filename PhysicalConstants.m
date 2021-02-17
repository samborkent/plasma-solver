classdef PhysicalConstants
    %PHYSICALCONSTANTS Class containing physical constants
    
    properties (Constant)
        BOLTZMANN = 1.38064852*10^-23;  % Boltzman constant [m2 kg / s2 K]
        AMU = 1.66053906660E-27;        % Atomic mass unit [kg]
        EV  = 1.602176634E-19;          % Electron volt [J]
        PLANCK = 6.62607004E-34;        % Planck constant [m2 kg / s]
        SOL = 299792458;                % Speed of light [m / s]
        
        % Atomic masses [kg]
        MASS_Sr = 87.62 * AMU;     % Strontium
        MASS_Li = 6.94 * AMU;      % Lithium
        MASS_Ti = 47.867 * AMU;    % Titanium
        MASS_O = 15.999 * AMU;     % Oxygen
    end
end

