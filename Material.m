classdef Material
    %MATERIAL Class holding material properties
    
    properties
        FORMULA
        ELEMENTS
        ENERGY_FORMATION
        ENERGY_BANDGAP
        LATTICE_A
        LATTICE_B
        LATTICE_C
        VOLUME
        DENSITY
    end
    
    methods
        function material = Material( formula, elements, E_f, E_bg, ...
                                        a, b, c, density )
            %MATERIAL Constructs and instance of this class
            material.FORMULA            = formula;
            material.ELEMENTS           = elements;
            material.ENERGY_FORMATION   = E_f;
            material.ENERGY_BANDGAP     = E_bg;
            material.LATTICE_A          = a;
            material.LATTICE_B          = b;
            material.LATTICE_C          = c;
            material.VOLUME             = a * b * c;
            material.DENSITY            = density;
        end
    end
end