classdef Material
    %MATERIAL Class holding material properties
    
    properties
        FORMULA
        ELEMENTS
        AMOUNT
        VOLUME
        DENSITY
        ENERGY_FORMATION
        ENERGY_TOTAL
    end
    
    methods
        function material = Material( formula, elements, amount, ...
                                      volume, density, E_f, E_tot )
            %MATERIAL Constructs and instance of this class
            material.FORMULA            = formula;      % Unit cell formula
            material.ELEMENTS           = elements;     % Array of elements in unit cell
            material.AMOUNT             = amount;       % Array of amount of each element in unit cell
            material.VOLUME             = volume;       % Unit cell volume [m^3]
            material.DENSITY            = density;      % Single crystal density [kg / m^3]
            material.ENERGY_FORMATION   = E_f;          % Formation energy per atom [J]
            material.ENERGY_TOTAL       = E_tot;        % Total energy of unit cell [J]
        end
    end
end