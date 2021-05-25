classdef PeriodicTable
    %PERIODICTABLE Class holding atomic properties
    % Source: PTable.com
    
    properties (Constant)
        % Lithium
        Li  = Atom( 3, ...                                  % Atomic number
                    'Li', ...                               % Symbol
                    'Lithium', ...                          % Name
                    6.94    * PhysicalConstants.AMU, ...    % Mass [kg]
                    167E-12, ...                            % Radius [m]
                    5.391   * PhysicalConstants.EV ...      % 1st ionization energy [J]
                   );
        % Oxygen
        O  = Atom( 8, ...                                   % Atomic number
                    'O', ...                                % Symbol
                    'Oxygen', ...                           % Name
                    15.999  * PhysicalConstants.AMU, ...    % Mass [kg]
                    48E-12, ...                             % Radius [m]
                    13.618  * PhysicalConstants.EV ...      % 1st ionization energy [J]
                   );
        % Argon
        Ar  = Atom( 18, ...                                 % Atomic number
                    'Ar', ...                               % Symbol
                    'Argon', ...                            % Name
                    39.948  * PhysicalConstants.AMU, ...    % Mass [kg]
                    71E-12, ...                             % Radius [m]
                    15.760  * PhysicalConstants.EV ...      % 1st ionization energy [J]
                   );
        % Titanium
        Ti  = Atom( 22, ...                                 % Atomic number
                    'Ti', ...                               % Symbol
                    'Titanium', ...                         % Name
                    47.867  * PhysicalConstants.AMU, ...    % Mass [kg]
                    176E-12, ...                            % Radius [m]
                    6.8280  * PhysicalConstants.EV ...      % 1st ionization energy [J]
                   );
        % Manganese
        Mn  = Atom( 25, ...                                 % Atomic number
                    'Mn', ...                               % Symbol
                    'Manganese', ...                        % Name
                    54.938  * PhysicalConstants.AMU, ...    % Mass [kg]
                    161E-12, ...                            % Radius [m]
                    7.4343  * PhysicalConstants.EV ...      % 1st ionization energy [J]
                   );
        % Strontium
        Sr  = Atom( 38, ...                                 % Atomic number
                    'Sr', ...                               % Symbol
                    'Strontium', ...                        % Name
                    87.62   * PhysicalConstants.AMU, ...    % Mass [kg]
                    219E-12, ...                            % Radius [m]
                    5.6952  * PhysicalConstants.EV ...      % 1st ionization energy [J]
                   );
        % Lanthanum
        La  = Atom( 57, ...                                 % Atomic number
                    'La', ...                               % Symbol
                    'Lanthanum', ...                        % Name
                    138.91  * PhysicalConstants.AMU, ...    % Mass [kg]
                    195E-12, ...                            % Radius [m]
                    5.577  * PhysicalConstants.EV ...       % 1st ionization energy [J]
                   );
    end
end

