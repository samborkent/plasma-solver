# Plasma Solver v5.0.0
Author: Tom Wijnands

Modified for two metal oxides by: Bart Boonstra

Reworked and modified for lithium metal oxides by: Sam Borkent

Based on 'Numerical modelling of the plasma plume propagation and oxidation
during pulsed laser deposition of complex oxide thin films', 2020, by
T. Wijnands, E.P. Houwman, G. Koster, G. Rijnders, and M. Huijben.

The original script by Tom Wijnands modeled the propagation of an PLD
plasma plume in 3D as result of ablation of a single crystal TiO2 target,
initial ablation is not included. Modified by Bart Boonstra to improve
performance and extended to include two metals (SrTiO3). Refactored and
extended by Sam Borkent to improve readability and performance and to
include two metals with highly different masses (Li4Ti5O12).

# Assumptions
* The target material is completely oxidized. Most oxide targets are 
    specified with a variable amount of oxygen in the unit cell (e.g.
    SrTiOx). This model assumes no oxygen vacancies or excess in the
    target. (e.g. SrTiO3, Li4Ti5O12)
* The excitation energy of the atoms in the unit cell are neglected
* Identical particles are assumed to seperate in space due to their
    initial velocity distribution, thus collisions between identical
    particles are neglected

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
* Field.veloBins    (1)
* Field.nParticles  (2)

# Change log
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

# Issues
* Number of particles in initialVelocityDistribution is
    a factor 10 too high

# To do
* Make fill inital plume and background matrix a function
* Implement folder selection into initialVelocityDistribution save feature
* Get correct formation energy of cubic spinel Li4Ti5O12 target
* Check which Li, Ti, and O molecules forms in the PLD plasma plume
* Check initial particle velocity distribution width
* Check restriction values
* Check ablation depth in target
* Check absorption coefficient value (Temp and photon energy dependend)
* Check crystal binding energy
* Check activation energies
* Check oxidation energies
* Check angular fitting parameter
* Check which data is relevant to print in config file

# Questions
* Why are all the array's in the Makestruct function double?
    Example: V = [0 20000 40000 0 20000 40000]
* Why is the crystal binding energy 5x the formation energy?
* What is the target density, and why is it 1?
* Why are O2 and Ti2O3 not included in mass array?
* Why is only the activation energy of TiO listed, and not of SrO or
  O2?
* Where does the 0.9 come from in Tom's average velocity?
* Why is the total number of ablated unit cells divided by the number of
    atoms in the unit cell?
