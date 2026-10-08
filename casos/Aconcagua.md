AoA: Angle of attack
Cd: Drag coefficient
Cl: Lift coefficient
Cm: Moment coefficient
# Purpose of simulation
Validate with the use of two CFD libraries (OpenFOAM and Nektar++) against  flight telemetry of ITBA Rocketry's  Aconcagua supersonic rocket. Libraries results will be compared against each other in order to evaluate which one gets more accuracy with same computing demands
Variables to be compared between simulations and telemetry are:
- Cd
- Cl
- Cm
Between simulations:
- DoF/error ratio for each variable
- Wall clock time
# Simulation structure 
Five base cases, each with the same AoA and progressively higher velocities were simulated in order to get as a result the variables curves, with the goal of comparing them against telemetry curves.
Velocities taken for the cases were 0.3, 0.6, 0.8, 1.2 and 1.8 Ma. 
## Geometry approximation and preparation

Holes, rail buttons, and engien

## Flow conditions
Aligning with thesis purpose of the reinforcement, verification and validation of a compressible flow solver, the compressible flow hypothesis was taken, as from 0.3 Ma onwards signs of variations in the density field of the fluid appear.
In consequence:
- 
# Solvers
## OpenFOAM: rhoCentralFoam

### Algorithm

### Associated problems

## Nektar++: CompressibleFlowSolver

Theory: 10.5 and 10.6 Karniadakis and Sherwin

It is well known the spectral methods degradation of convergence and monotonicity due to Gibbs phenomenon in transonic and supersonic applications.

### Algorithm


# Preprocessing
## Rocket model
The .stl file of the Aconcagua is extracted from the full model made in Fusion. Geometry is extracted with the highest resolution in order to get best result possible when meshing.
## Mesh definition
### Symmetry simplifications (TBR)
Full symmetry was used in order to capture compressible related phenomena properly.
### Inlet and outlet
