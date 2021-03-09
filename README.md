# Plasma Solver v5.0.0
Author: Tom Wijnands

Rewritten and modified for targets of any composition: Sam Borkent

Based on 'Numerical modelling of the plasma plume propagation and oxidation
during pulsed laser deposition of complex oxide thin films', 2020, by
T. Wijnands, E.P. Houwman, G. Koster, G. Rijnders, and M. Huijben.

The original script by Tom Wijnands modeled the propagation of an PLD
plasma plume in 3D as result of ablation of a single crystal TiO2 target,
initial ablation is not included. Rewritten and extended by Sam Borkent
to improve readability and performance and to support targets of any
composition.

# Assumptions
* The excitation energy of the atoms in the unit cell are neglected.
* Identical particles are assumed to seperate in space due to their
    initial velocity distribution, thus collisions between identical
    particles are neglected.
* The background gas starts out with zero velocity, as the room tempreature
    velocity is aroung 400 m/s, which is two orders smaller than the plume
    particle velocities. The particles in the background gas don't have a
    preferential direction, so their average velocity zeros out.
* Only head-on collisions are included, so there is no net exchange of
    particles between angular bins.
* Particles can only gain positive momentum, when a lighter plasma particle
    collides with a heavier background gas particle its momentum is set to
    zero, and the background gas particle gains a little momentum, so
    momentum is not conserved.
* Collisions are only possible if the background particles move slower than
    the plume particles

# Naming conventions
Loosly resemble paper names and follow MATLAB style guidelines v1.3

* radius    : Radius / Distance
* velo      : Velocity
* dist      : Distribution
* init      : Initial
* n         : Number of particles
* bg        : Background gas
* uc        : Unit cell
* bin       : Computational bin with similar position, velocity, etc.

Field names and corresponding index value: (see Field enumeration class)
* Field.nParticles  (1)     Number of particles
* Field.nCollisions (2)     Number of collisions
* Field.O           (1)     Oxygen

# Issues
* When making crazy composite targets, lithium gain too much energy

# Change log

09-03-21:
* Worked on main loop, made start with collisions, no propagation yet

06-03-21:
* Restructure main loops to only excecute code when necessary. Angles get
    skipped if number of particles per angle is smaller than minimum
    number of particles per angle.

04-03-21:
* Implemented target density
* Implemented mixture targets

27-02-21:
* Added memory check to prevent program from running if not enough memory
    is available

25-02-21:
* Changed initialVelocityDistribution to output all particle velocity
    distributions in one matrix
* Changed Fields, matrices do not have to hold velocity axis, only number
    of particles in each velocity bin, and the average number of collisions
* Added dimension to plume matrix to hold all species simultaneously

24-02-21:
* Changed file structure
* Added documentation

23-02-21:
* Started work on the main program, replaced struct with a 3D matrix to
    improve performance
* Added Field enumerator class to index fields withing the 3D matrix to
    improve readability

22-02-21:
* Added amount array and number of atoms property to Material, replaced
    elements string for PeriodicTable reference
* Fixed issue where velocity changes with number of atoms. Number of atoms
    still factor 10 to high

20-02-21:
* Edited initialVelocityDistribution to try to fix of shift with increaed
    atoms
* Made Material class, and UnitCells class to hold all material properties

19-02-21:
* Added Atom class holding Atom properties, added PeriodicTable class
    holding all atoms
* Improved kinetic energy term in initialVelocityDistribution and added
    excitation energy
* Multiplied nParticleAngle by number of atoms in unit cell

18-02-21:
* Added added atomic radii and Mn to PhysicalConstants
* Added laser fluence variable
* Fixed initial velocity distribution to reproduce results of paper
* Added notes where the there is a discrepancy between paper and script

17-02-21:
* Edited to include Git support
* Changed naming conventions to more closely resemble the names used in
  the paper and to follow MATLAB style guidelines
* Changed maximum velocity to include limits of Li
* Wrote the initial velocity distibution as a function
* Wrote the configuration file creator as a function
* Changed the mass array for LTO and included complex compounds
* Added relavant material properties to PhysicalConstants class
* Move unused code to oldCode

16-02-21:
* Work on initial velocity distribution

12-02-21:
* Started refactoring main loop
* Replaced looping through angles to index looping
* Refactored collision cross section

08-02-21:
* Starting top-to-bottom to refactor the code to increase readability
* Adjust naming conventions and replace Dutch variable names
* Change 1/0 bools into true/false bools
* Denoting lines with unclear purpose with %?%
* Reorder variables and place them in more clear groups
* Added distiction between lower/Camel/UPPER case for
  variables/Parameter/CONSTANTS
* Made make dir conditional so it is skipped if dir already exists
* Add date and time to config file

# To do
* Make the model work without collisions or oxidation
* Replace divisions that are performed often with multiplications
* Remove code that is unnecessary
* Implement excitation energy
* Preacclocate object arrays
* Add correcting to unit cell volume for non-right angles
* Make fill inital plume and background matrix a function
* Implement folder selection into initialVelocityDistribution save feature
* Get correct formation energy of cubic spinel Li4Ti5O12 target
* Get materials data straight from PeriodicTable and MaterialsProject
    websites
* Atomatically determine non-ionized oxidation states in plasma plume

# To validate
* Species that form in plasma plume
* Angular fitting parameter
* Initial particle velocity distribution width
* Restriction values
* Ablation depth in target
* Absorption coefficient (Temp and photon energy dependend)
* Activation energies
* Oxidation energies
* Which data to print in config file

# Questions
* Why does the first angle bin not have the largest number of particles?
* Why are all the array's in the Makestruct function doubled?
    Example: V = [0 20000 40000 0 20000 40000]
* Why are O2 and Ti2O3 not included in mass array?
* Why is only the activation energy of TiO listed, and not of SrO or
  O2?
