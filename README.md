# Lane-Free Dynamic Programming / DDP Path Planning

This repository contains research code and supporting material for path-planning of automated vehicles in lane-free traffic. The code is associated with constrained Dynamic Programming / Differential Dynamic Programming (DDP) formulations for computing safe, smooth, and efficient trajectories in the presence of road boundaries and obstacle vehicles.

The repository is related to the following work:

> Panagiotis Typaldos and Markos Papageorgiou,  
> **A Constrained Differential Dynamic Programming Algorithm for Automated Vehicles in Lane-Free Traffic**,  
> MFTS 2024.

## Overview

The path-planning problem is formulated as an optimal control problem for an automated vehicle moving in a lane-free road environment. The objective function balances:

- vehicle advancement toward a desired speed,
- passenger comfort through acceleration penalties,
- obstacle avoidance,
- road-boundary satisfaction,
- smooth longitudinal and lateral motion.

The vehicle dynamics are modeled with double-integrator equations in the longitudinal and lateral directions. Lateral acceleration bounds are used to keep the vehicle within road boundaries, while obstacle-avoidance terms guide the vehicle around surrounding vehicles.

## Status

This repository contains research code associated with published and in-progress work. It is maintained as a reference implementation and experiment archive rather than a production software package.

## Repository map

```text
mathematica/
  Supporting derivations, symbolic calculations, and exploratory notebooks.

project MFTS2024/
  Main code and experiment files associated with the MFTS 2024 constrained DDP work.

sim/
  Simulation and visualization utilities for lane-free path-planning experiments.

outputs/
  Generated outputs from local runs. These files are treated as experiment artifacts.