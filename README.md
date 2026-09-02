# BESFEM: Battery Electrode Simulation using MFEM

BESFEM (**B**attery **E**lectrode **S**imulation using **MFEM**) is a high-performance finite element framework for simulating lithium-ion battery electrodes. Built on top of **MFEM**, **MPI**, and **HYPRE**, BESFEM enables parallel electrochemical simulations of realistic battery microstructures using the **Smoothed Boundary Method (SBM)**.

The framework supports both **half-cell** simulations, multiple active material chemistries, and particle-resolved modeling using either diffusion or Cahn–Hilliard-based transport models. The **full-cell** model is currently under construction. 

---

# Project Structure

```text
BESFEM/
│
├── include/             # Header files
├── src/                 # Source files
├── inputs/
│   ├── mesh/            # Mesh and TIFF geometries
│   ├── materials/       # Material property tables
│   └── run_config.txt   # Simulation configuration
│
├── outputs/             # Simulation outputs
├── plotting/            # Visualization notebooks
├── tests/               # Unit tests
├── docs/                # Generated Doxygen documentation
└── bin/                 # Compiled executable
```

---

# Building BESFEM

Ensure that MFEM, MPI, and HYPRE are installed and available.

```bash
# Clone the repository
git clone https://github.com/HCY-Group/BESFEM.git

# Enter the project
cd BESFEM

# Load required modules (HPCC example)
module load MFEM

# Compile
make

# Executable
cd bin
```

> Depending on your system, you may need to update the makefile with the appropriate MFEM and HYPRE include/library paths.

---

# Running BESFEM

All simulations are configured through

```text
inputs/run_config.txt
```

Run using the default configuration:

```bash
mpirun -np 8 ./battery_simulation
```

Specify a configuration file:

```bash
mpirun -np 8 ./battery_simulation -cfg ../inputs/run_config.txt
```

Save terminal output:

```bash
mpirun -np 8 ./battery_simulation -cfg ../inputs/run_config.txt > output.txt
```

Extract timestep information:

```bash
grep timestep output.txt > timestep.txt
```

---

# Simulation Workflow

```text
Geometry (TIFF to Mesh)
          │
          ▼
Define Domain Parameters (SBM)
          │
          ▼
Assign Materials
          │
          ▼
Initialize Concentrations & Potentials
          │
          ▼
Solve Coupled Electrochemical Equations
          │
          ▼
Write Output Files
```

---

# Configuration

All simulation settings are specified in `inputs/run_config.txt`.

A typical configuration is shown below.

```ini
mode = half
electrode = cathode
mesh_file = ../inputs/mesh/colored_labels_labels.tif

combine_particles = true

row_begin = 0
row_end = 100
column_begin = 0
column_end = 100

amr_levels = 1
coarsen_factor = 2

dt = 1.2e-05
dh = 8.0e-07
gc = 6.38e-12

stop_mode = steps       
VCut = 3.41
num_steps = 1501

Cr = 1.0
Vsr0 = 0.9466

cathode_materials = LFP
init_cathode_particles = 0.043

init_CnE = 0.001

init_BvC = 3.359
init_BvE = -0.1
```

---

### Simulation

* `mode`

  * `half`
  * `full`

* `electrode`

  * `anode`
  * `cathode`

---

### Geometry

* `mesh_file`

  * TIFF geometry

* `combine_particles`

  * `true` — treat all particles as a single particle.
  * `false` — solve each particle independently.

* TIFF Crop Bounds

  * `row_begin`
  * `row_end` 
  * `column_begin` 
  * `column_end`

* `amr_levels`

  * AMR is supported for 1 level of refinement.

* `coarsen_factor`
  * Select a factor at which the grid will coarsen before AMR usage. Typically 2 or 4. 


---

### Stopping Criteria

* `stop_mode`
  * Choose how the simulation will stop - either by steps or voltage.

* `num_steps`

  * Total number of timesteps.

* `VCut`
  * Cut Off Voltage.

---

### Materials

Materials are assigned in the same order as the particle groups.

Example

```ini
cathode_materials = LFP,LFP,NMC
```

assigns

* Particle 0 → LFP
* Particle 1 → LFP
* Particle 2 → NMC

Currently supported materials include

**Cathodes**

* LFP
* NMC

**Anodes**

* Graphite
* Carbon

---

### Initial Conditions

Particle lithiation is specified per particle.

Example

```ini
init_cathode_particles = 0.15,0.20,0.10
```

Electrolyte concentration

```ini
init_CnE = 0.001
```

Initial electrode and electrolyte potentials

```ini
init_BvA = -0.10
init_BvC = 3.40
init_BvE = -0.10
```

The battery discharge and charge rate can be declared. 
A positive rate (+) is a discharge simulation and a negative rate (-) is a charge simulation.

```ini
Cr = 1.0     # discharge simulation
Cr = -1.0    # charge simulation
```

---

### Output

Simulation results are written to the `outputs/` directory and include quantities such as

* Particle concentration
* Electrolyte concentration
* Electrode potentials
* Electrolyte potential
* SBM phase fields
* Mesh files
* Simulation logs

These outputs may be visualized using **PyGLVis** or other MFEM-compatible visualization tools.

---
### Boundary Conditions

BESFEM uses a fixed orientation convention for the external boundaries of the computational domain. Boundary attribute indicators are assigned as follows:

| Boundary | Indicator |
| -------- | --------: |
| South    |       `0` |
| East     |       `1` |
| North    |       `2` |
| West     |       `3` |

The geometry is oriented so that the **anode current collector is located on the west side** of the domain and the **cathode current collector is located on the east side**. The boundary conditions are assigned in `src/BoundaryConditions.cpp`. 

For a full-cell simulation, the electrode arrangement is therefore

```text
West                                                    East
Anode current collector → Anode → Electrolyte → Cathode → Cathode current collector
```

The solid-phase electrical-potential boundary conditions are applied at the corresponding current collectors:

* The anode solid potential is constrained on the west boundary, with boundary indicator `3`.
* The cathode solid potential is constrained on the east boundary, with boundary indicator `1`.

The north and south boundaries, with indicators `2` and `0`, respectively, represent the transverse edges of the geometry. Unless otherwise specified by a particular model, these boundaries use no-flux or electrically insulating boundary conditions.

For half-cell simulations, the current collector remains on the electrode-specific side:

* Anode half-cell: current collector on the west boundary.
* Cathode half-cell: current collector on the east boundary.

The electrolyte-facing boundary is located on the side opposite the current collector. Internal active-material/electrolyte interfaces are represented through the SBM phase fields and are not assigned external mesh-boundary indicators. Electrochemical reaction terms are evaluated along these diffuse internal interfaces.

In full-cell geometries, the electrolyte phase is retained only when it forms a connected pathway that contacts both the anode and cathode phases. Isolated electrolyte regions that do not connect both electrodes are removed from the electrolyte connectivity mask.


---

# Generating Documentation

Generate the Doxygen documentation:

```bash
module load Doxygen
module load Graphviz

doxygen Doxyfile
```

Open the generated documentation

```bash
cd docs/html
```

Preview locally

```bash
python3 -m http.server 8001 --bind 127.0.0.1
```

Then open

```
http://127.0.0.1:8001
```

---

# Visualization

Install PyGLVis

```bash
pip install glvis
```

Example visualization notebooks are located in

```text
plotting/
```

Update the mesh and GridFunction filenames within the notebook to visualize different simulation outputs.

---

# Governing Equations

## Cathode Concentration

```math
\frac{\partial C_c}{\partial t}
=
\frac{1}{\psi_c}
\nabla\cdot
\left(
\psi_c D_c \nabla C_c
\right)
-
\frac{|\nabla\psi_c|}{\psi_c}r_c
```

---

## Cathode Potential

```math
\nabla\cdot
\left(
\psi_c\kappa_c\nabla\phi_c
\right)
-
|\nabla\psi_c|z_-Fr_c
=
0
```

---

## Electrolyte Concentration

```math
\frac{\partial C_e}{\partial t}
=
\frac{1}{\psi_e}
\nabla\cdot
\left(
\psi_eD_e\nabla C_e
\right)
+
\frac{|\nabla\psi_c|}{\psi_e}
\frac{r_ct_-}{\nu_+}
+
\frac{|\nabla\psi_a|}{\psi_e}
\frac{r_at_-}{\nu_+}
```

---

## Electrolyte Potential

```math
\nabla\cdot
\left[
\psi_e
(z_+m_+-z_-m_-)
FC_e
\nabla\phi_e
\right]
+
|\nabla\psi_c|
\frac{r_c}{\nu_+}
+
|\nabla\psi_a|
\frac{r_a}{\nu_+}
=
\nabla\cdot
\left[
\psi_e
(D_- - D_+)
\nabla C_e
\right]
```

---

## Anode Concentration

```math
\frac{\partial C_a}{\partial t}
=
\frac{1}{\psi_a}
\nabla\cdot
\left[
\psi_a
M_a
\nabla
\left(
\frac{\partial f_G}{\partial C_a}
-
\varepsilon\nabla^2C_a
\right)
\right]
-
\frac{|\nabla\psi_a|}{\psi_a}r_a
```

---

## Anode Potential

```math
\nabla\cdot
\left(
\psi_a\kappa_a\nabla\phi_a
\right)
-
|\nabla\psi_a|z_-Fr_a
=
0
```

---

# Citation

If you use BESFEM in your research, please cite the appropriate publication describing the framework.

---

# License

This project is released under the license included with the repository.
