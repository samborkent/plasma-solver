classdef Material
    %MATERIAL Class holding material properties
    
    properties
        FORMULA
        ELEMENTS
        AMOUNT
        N_ATOM
        ENERGY_FORMATION
        LATTICE_A
        LATTICE_B
        LATTICE_C
        VOLUME
        DENSITY
        ATOM_DENSITY
    end
    
    methods
        function material = Material( formula, elements, amount, ...
                                        E_f, a, b, c, density )
            %MATERIAL Constructs and instance of this class
            material.FORMULA            = formula;
            material.ELEMENTS           = elements;
            material.AMOUNT             = amount;
            material.N_ATOM             = sum(amount);
            material.ENERGY_FORMATION   = E_f;
            material.LATTICE_A          = a;
            material.LATTICE_B          = b;
            material.LATTICE_C          = c;
            material.VOLUME             = a * b * c;
            material.DENSITY            = density;
            material.ATOM_DENSITY       = material.N_ATOM / material.VOLUME;
        end
    end
end