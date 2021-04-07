# Plasma Solver v5.0.0
Author: Tom Wijnands

Rewritten and modified for targets of any composition: Sam Borkent

Based on 'Numerical modelling of the plasma plume propagation and oxidation
during pulsed laser deposition of complex oxide thin films', 2020, by
T. Wijnands, E.P. Houwman, G. Koster, G. Rijnders, and M. Huijben.

The original script by Tom Wijnands modeled the propagation of an PLD
plasma plume in 2D as result of ablation of a single crystal TiO2 target,
initial ablation is not included. Rewritten and extended by Sam Borkent
to improve usability and performance and to support targets of any
composition consisting of species of any mass.

# Assumptions
* The excitation energy of the atoms in the unit cell induced by the
	ablation laser is neglected.
* The background gas starts out with zero velocity: the room temperature
    velocity is aroung 400 m/s, which is two orders smaller than the plasma
    particle velocities (~10^4 m/s). The particles in the background gas do
    not have a preferential direction, so their net velocity is zero.
* All collisions are head-on: collisions at random angles in all directions
    cancel out.
* There is no net exchange of particles between angular bins: the number of
    particles entering and leaving each angular bin each time step is
    similar.
* Particles can only gain positive momentum. When a light particle collides
    with a heavier particle, its momentum is set to zero, and the heavy
    particle gains a little momentum, so momentum is not conserved.
* Collisions of species with the same mass are neglected, as for
    completely elastic collisions they would interchange velocities,
    resulting in no net change of density.

# Naming conventions
Follow MATLAB style guidelines v1.3

* bg      : Background
* bin     : Computational bin with similar position and velocity
* n       : Number of particles/bins
* velo    : Velocity
* uc      : Unit cell

# Issues
* All background particles in the first few bins get scattered due to a
    much higher density of plasma particles (10^15) in the first radial bin
    compared to background particle density (10^12). So, even with a
    low collision rate (e.g. 0.001), all background particles get scattered,
    resulting in vacuum formation in front of the target. The number of
    particles seems to be conserved, even though the area under the density
    curve of the propagating background gas is not conserved.
* When making complex composite targets, lithium gains too much energy.
* Cannot reproduce Tom's results: particles travel faster right after
    after ablation but propagate slower.

# Questions
* Why does the first angle bin not have the largest number of particles?

# To Do

Main functionality
* Implement proper smoothing/fitting of data
* Verify results with one heavy metal collision
* Counting the number of collisions per bin
* Reproduce the Ti 1D propagation plot of the paper
* 2D density plot
* Collisions for two species
* Collisions for light metals (Li)
* Collisions for any number of species with any mass
* Oxidation of plasma species
* Collisions between unoxidized and oxidized species

Physical improvements
* Implement excitation energy in initialVelocityDistribution
* Calculate heat dissipation into target
* Include a temperature gradient origination from the heated substrate,
    resulting in a lower density and higher kinetic energy of background
    species near the substrate.

Additional features
* Add correction to unit cell volume for non-right angles
* Implement folder selection into initialVelocityDistribution save feature
* Get materials data straight from P-Table and MaterialsProject websites
* Atomatically determine non-ionized oxidation states in plasma plume

Clean-up
* initialVelocityDistribution
* Rename nPlasmaParticlesPerRadius to nParticlesPerRadius
* Include plotting of 1D propagation in nParticlesPerRadius function

Performance
* Replace all arrays with values in range 10^(+-38) with single-precision
* Replace find with logical indexing where possible
* Replace division with multiplication where possible
* Include second pass to update background particles in the main loop

Validation
* Which species form in plasma
* Angular fitting parameter
* Initial particle velocity distribution width
* Ablation depth in target
* Absorption coefficient (Temp and photon energy dependend)
* Activation energies
* Oxidation energies
* Which data to print in config file
* Formation energy of cubic spinel Li4Ti5O12 target

# Change log

06-04-21:
* Cleaned up readme
* Worked on number of collisions per bin. It's working, but results are
    wrong.
* Tried to implement smoothing/fitting.

26-03-21:
* Made a second pass per time step to update background particle positions,
    to ensure all particles get updated.
* Tried to fix collisions. It's still not working as planned. All bg
    particles close to target get scattered due to the high initial plasma
    particle density, creating a vacuum near the target.
* Worked on number of collision system, number of collisions does not
    update properly.

25-03-21:
* Placed the particle update after the collision calculation. Still need to
    have the distance traveled be different for collided particles and to
    include collisions in the starting bin.
* Tried to implement number of collisions, started a new main file.
    (previous code in testPlasmaSolver, even older code in
    testPlasmaSolverOld) Added another dimension to the matrix to count
    collisions. It's (kind of) working, but does not give expected results.

20-03-21:
* The velocities of all particles in the last radial bin gets set to zero.
* Fixed updating of matrices.

18-03-21:
* Implemented first working version of collision with one species and
    background. Still needs to be validated.
* Added debug messages (labeled [DEBUG]) for mean collision probability,
    collision probability per bin, and number of particles.
* Included normalization facotr for velocity weight factors.
* Added some default values.

12-03-21:
* Made a collisionProbability function to try to make things work.

10-03-21:
* Renamed the current .m file to testPlasmaSolver.m, and created 
    plasmaSolver.m which will only hold my own working code.
* Managed to kind of reproduce 1D propagation, although with quantization
    problem. And particles don't all go to the last radial bin eventually.
* Created function to initialize and fill background matrix (fillBGMatrix).
* Created nPlasmaParticlesPerRadius function that sums particles of all
    velocity bins within current radial bin
* Added envelope for plotting 1D propagation, which is normalized to
    conserve the number of particles
* Created angularDistribution function to calculate angular distribution
    of plasma particles
* Did rudementary collisions. It works, but the number of collisions is
    too high and particles come to a stand still. Also background density
    is assumed constant, background gas is not modeled yet.

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
