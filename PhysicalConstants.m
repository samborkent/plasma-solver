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
    end
end
