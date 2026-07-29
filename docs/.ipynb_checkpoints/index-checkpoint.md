---
title: Home
layout: home
nav_order: 1
---

## Battery Electrode Simulation using MFEM

Battery Electrode Simulation Finite Element Method (BESFEM) is a high-performance finite element framework for simulating electrochemical processes in lithium-ion battery electrodes.

It is built on MFEM, MPI, and HYPRE so that it supports parallel simulations of battery microstructures.

---

## About BESFEM

The performance of battery electrode is highly dependent by the geometry, material composition and behavior of individual particles within the electrode microstructure. Electrode microstructures can be provided in the form of TIFF images, and each individual particle clusters are identified with distinct materials and initial conditions. BESFEM enables the investigation of such effects using particle-based finite element analysis. The framework models the computational analysis of electrochemical processes in batteries.

BESFEM currently supports half-cell simulations for anode and cathode systems. Development of the full-cell model is ongoing.

---

## Key Capabilities

### Microstructure-resolved simulation

BESFEM characterizes microstructures of battery electrodes and retains data regarding individual particle clusters.

### Smoothed Boundary Method

The Smoothed Boundary Method is used to characterize material interface boundaries inside the finite element model.

### Particle-level material assignment

Different active materials and initial lithiation values can be assigned to individual particle groups within the input microstructure.

### Multiple transport models

Particle transport can be modeled using diffusion-based or Cahn–Hilliard-based formulations.

### Parallel finite element computing

BESFEM is built with MFEM, MPI, and HYPRE to support parallel numerical simulations.

### Adaptive mesh refinement

The system supports adaptive mesh refinement to enhance the resolution of specific regions in the microstructure.

---

## Supported Systems

BESFEM currently supports half-cell simulations for both cathode and anode electrodes.

Supported active materials include:

| Electrode | Materials                                                            |
| :-------- | :------------------------------------------------------------------- |
| Cathode   | Lithium iron phosphate (LFP) and nickel manganese cobalt oxide (NMC) |
| Anode     | Graphite                                                             |


The materials are allocated based on the particle groups defined in the microstructure input.

---

## Simulation Workflow

A typical BESFEM simulation follows these stages:

1. Provide a TIFF image describing the electrode microstructure.
2. Generate the computational mesh and Smoothed Boundary Method phase fields.
3. identify the individual particle groups in the microstructure.
4. Assign active materials and initial conditions.
5. Configure the simulation through `run_config.txt`.
6. Solve the coupled electrochemical equations.
7. Write the concentration, potential, mesh, and simulation output files.
8. Visualize and analyze the results.

---

## Simulation Configuration

BESFEM simulations are controlled through a configuration file. Users can specify settings such as:

* the electrode being simulated
* the microstructure input file
* active materials
* particle-level initial lithiation
* electrolyte concentration
* initial electrode and electrolyte potentials
* simulation stopping criteria
* the number of timesteps
* voltage cutoff
* particle grouping
* adaptive mesh refinement levels.

This configuration-based workflow allows simulations to be changed without modifying the source code.

---

## Simulation Outputs

BESFEM produces results describing the spatial and temporal behavior of the simulated electrode, including:

* active-particle concentration
* electrolyte concentration
* electrode potential
* electrolyte potential
* Smoothed Boundary Method phase fields
* finite element mesh files
* simulation logs.

The generated finite element results can be visualized using PyGLVis and other MFEM compatible visualization tools.

---

## Documentation

Use the navigation menu to access the BESFEM API documentation and other available project resources.

The API reference is generated from the source code using Doxygen and provides detailed information about BESFEM classes, functions, files, and data structures.

---

## Project Status

BESFEM is under active development.

The current implementation supports half-cell electrode simulations. The full-cell model and additional simulation capabilities are still being developed.

---

## Contributing

Contributions that improve the source code, testing, documentation, numerical methods, or simulation capabilities are welcome.

Please use the BESFEM GitHub repository to report issues, propose changes, or submit a pull request.

---

## Citation

When using BESFEM in published research, please cite the appropriate publication or software release describing the framework.
