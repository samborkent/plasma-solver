classdef Atom
    %ATOM Class holding atomic properties
    
    properties
        NUMBER
        SYMBOL
        NAME
        MASS
        RADIUS
        ENERGY_FI
    end
    
    methods
        function atom = Atom(number, symbol, name, mass, radius, energy_fi)
            %ATOM Construct and instance of this class
            atom.NUMBER     = number;
            atom.SYMBOL     = symbol;
            atom.NAME       = name;
            atom.MASS       = mass;
            atom.RADIUS     = radius;
            atom.ENERGY_FI  = energy_fi;    % First ionization energy [J]
        end
    end
end