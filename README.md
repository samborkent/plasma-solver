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
include two metals with highly different masses (Li2TiO3).

# Change log:

17-02-21:
  * Edited to include Git support
  * Changed naming conventions to more closely resemble the names used in
      the paper and to follow MATLAB style guidelines

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
  * Check with data is relevant to print in config file
  * Improve naming conventions following paper and Matlab guidelines
  * Correct inital velocity distribution ( f_i(v) )
  * Implement excitation energy

# Questions:

  * Why is the crystal binding energy 5x the formation energy?
  * It says that is 'Oxidation' is set to 1 (true) there is no oxidation
      possible, is this correct?
  * What is the target density, and why is it 1?
  * Why are O2 and Ti2O3 not included in mass array?
  * Why is only the activation energy of TiO listed, and not of SrO or
      O2?
  * What do the three constants in the initial plume settings block do?
  * Why clear variables from memory within loop, instead of overwriting?
  * Where does the 0.9 come from in Tom's average velocity?
  * Why is the average velocity twice divided by the number of atoms in
      the unit cell?
  * V_Delta is experimentally fitted according to the paper, but is
      exactly 5?

# Naming conventions:
Loosly resemble paper names and follow MATLAB style guidelines v1.3

  * lowers case : variable  (can be altered within script)
  * Camel Case  : Parameter (can be altered by user before running script)
  * UPPER CASE  : CONSTANT  (can not be altered)
  * t         : Time
  * x         : Distance
  * rad       : Angle
  * n         : Number of particles
  * v         : Velocity
  * k         : Number of collisions
  * p         : Pressure
  * temp      : Temperature
  * e         : Energy
  * mo        : Moving particles
  * st        : Static particles
  * bg        : Background gas particles
  * col       : Colliding particles
  * ncol      : NOT colliding particles 
  * m         : Current space matrix of all plume particles
  * m_1       : New space matrix of all plume particles
  * m_bg      : Current space matrix of all background gas particles
