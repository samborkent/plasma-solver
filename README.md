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

# Change log:

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

# To Do:
* Fix problem where velocity distribution shifts to the right if the number
    of atoms in the unit cell are increased (probably velocityAverage)
* Get correct lattice parameters, unit cell volume, and formation energy
    of spinel Li4Ti5O12 target
* Check which Li, Ti, and O molecules forms in the PLD plasma plume
* Check initial particle velocity distribution width
* Check restriction values
* Figure out the purpose of lines marked with %?%
* Check ablation depth in target
* Check adsorption coefficient value (Temp and photon energy dependend)
* Check crystal binding energy
* Check activation energies
* Check oxidation energies
* Check angular fitting parameter
* Check inital maximal velocity
* Check velocity delta
* Check which data is relevant to print in config file
* Implement excitation energy

# Questions:
* Why is the crystal binding energy 5x the formation energy?
* What is the target density, and why is it 1?
* Why are O2 and Ti2O3 not included in mass array?
* Why is only the activation energy of TiO listed, and not of SrO or
  O2?
* Why clear variables from memory within loop, instead of overwriting?
* Where does the 0.9 come from in Tom's average velocity?
* Why is the average velocity twice divided by the number of atoms in
  the unit cell?
* Why is the total number of ablated unit cells divided by the number of
    atoms in the unit cell?

# Naming conventions:
Loosly resemble paper names and follow MATLAB style guidelines v1.3

* radius    : Radius / Distance
* velo      : Velocity
* n         : Number of particles
* bg        : Background gas
* uc        : Unit cell

# Assumptions
* The target material is completely oxidized. Most oxide targets are 
    specified with a variable amount of oxygen in the unit cell (e.g.
    SrTiOx). This model assumes no oxygen vacancies or excess in the
    target. (e.g. SrTiO3, Li4Ti5O12)
