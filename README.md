# Plasma Solver
Author: Tom Wijnands, Sam Borkent

Based on 'Numerical modelling of the plasma plume propagation and oxidation
during pulsed laser deposition of complex oxide thin films', 2020, by
T. Wijnands, E.P. Houwman, G. Koster, G. Rijnders, and M. Huijben.

The original script by Tom Wijnands modeled the propagation of a PLD
plasma plume in 2D as result of ablation of a single crystal TiO2 target;
initial ablation is not included. Rewritten and extended by Sam Borkent
to improve usability and performance and to support targets of any
composition containing species of any mass.

# Assumptions

* The laser spot is assumed to be square, resulting in a circular plasma
    cross-section.
* The background gas starts out with zero velocity. The room temperature
    velocity is around 400 m/s, which is two orders smaller than the plasma
    particle velocities (~10^4 m/s). The particles in the background gas do
    not have a preferential direction, so their net velocity is zero.
* All collisions are head-on: collisions at random angles in all directions
    cancel.
* There is no net exchange of particles between angular bins: the number of
    particles entering and leaving each angular bin each time step is
    similar.
* Collisions of species with the same mass are neglected, as for
    completely elastic collisions they would interchange velocities,
    resulting in no net change of density.
* Particles are assumed static if their velocity is lower than the velocity
    required to travel one radial bin in one time step.
* Oxidation only occurs with the background gas, not with plasma
    particles. In non-oxygen background gasses the material oxidizes
    at the substrate surface solely from oxygen originating from the
    target.
* A sole particle can only undergo one collision per time step, which is valid for
    sufficient resolution if collisions between plasma particles are neglected.
    
# Changes compared to Wijnands's model

* Background gas particles can have possess every velocity and can be modeled
    similar to plasma species.
* Double matrices are used for all data instead of structs, which improves
    performance and makes it easier to acces specific data.
* The code was written to be general purpose, so it works for any material,
    instead of only for TiO2.
* Adjusted collision rate calculation, and included a velocity weight term,
    so collisions with slower moving particles are more likely than
    collisions with faster moving particles.
* Implemented multiple initial velocity distributions, including a Maxwell-
    Boltzmann distribution, which is more physically accurate than a Gaussian
    distribution.

# Naming conventions

Follow MATLAB style guidelines v1.3

* bg      : Background
* bin     : Computational bin with similar position and velocity
* n       : Number of particles/bins
* velo    : Velocity
* uc      : Unit cell
