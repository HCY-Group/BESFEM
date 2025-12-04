# BESFEM: Battery Electrode Simulation using MFEM

BESFEM is a high-performance, MPI-enabled electrochemical simulation framework for lithium-ion battery electrodes.  
It supports **half-cell** and **full-cell** simulations, **Cahn–Hilliard** and **diffusion-based** transport models, and the **Smoothed Boundary Method (SBM)** for diffuse interfaces.  
The code is written in C++ and built on top of the **MFEM** finite-element library.


---

## Project Structure

```
BESFEM/
│
├── include/             # Header files for physics modules
├── src/                 # Source files for all simulations
├── inputs/
│   ├── mesh/            # Mesh files
│   ├── distance/        # SBM distance fields
│   └── constants/       # Material + parameter files
│
├── outputs/
│   └── Results/         # Auto-generated simulation outputs
│
├── tests/               # Unit tests
├── plotting/            # Plotting files
└── bin/                 # Compiled executables
```

---

## Building BESFEM

Ensure MFEM, HYPRE, and MPI (OpenMPI or MPICH) are installed and available.

```bash
# clone the repository
git clone --single-branch --branch AB_OOP https://gitlab.msu.edu/hcy/besfem.git

# enter into BESFEM folder
cd BESFEM

# compile all of the code - you may need to update  the makefile MFEM/HYPRE include + library paths as needed
make 

# create directory structure for outputs
mkdir outputs
cd outputs
mkdir Results

# enter folder with executable file
cd bin
```

---

## Running Simulations

### Full Cell Example
```bash
mpirun -np 2 ./battery_simulation \
    -mode full \
    -m ../inputs/mesh/Mesh_40x60x3_3D_disk_full.mesh \
    -dA ../inputs/distance/dsF_A_40x60x3_3D_disk_full.txt \
    -dC ../inputs/distance/dsF_C_40x60x3_3D_disk_full.txt \
    -t ml
```

### Half Cell Example (Cathode)
```bash
mpirun -np 4 ./battery_simulation \
    -mode half \
    -elec cathode \
    -m ../inputs/mesh/Mesh_40x60_F00.mesh \
    -dC ../inputs/distance/dsFC_41x61_F00.txt \
    -t ml
```

### Half Cell Example (Anode)
```bash
mpirun -np 3 ./battery_simulation \
    -mode half \
    -elec anode \
    -m ../inputs/mesh/Mesh_40x60_F00.mesh \
    -dA ../inputs/distance/dsFA_41x61_F00.txt \
    -t ml
```

---

## Command Line Options

| Option                  | Description                             |
| ----------------------- | --------------------------------------- |
| `-m <MeshFile>`         | Path to `.mesh` file                    |
| `-mode <half/full>`     | Select simulation mode                  |
| `-elec <anode/cathode>` | Required for half-cell mode             |
| `-dA <file>`            | Anode distance field (`.txt`)           |
| `-dC <file>`            | Cathode distance field (`.txt`)         |
| `-o <order>`            | Finite element polynomial order         |
| `-t <ml/v>`             | Mesh type: MATLAB (`ml`) or voxel (`v`) |
| `-n <steps>`            | Number of time steps                    |

---

## Generating Doxygen Documentation

```bash
module load Doxygen
doxygen Doxyfile
cd html
```
Preview doxygen locally:
```bash
python3 -m http.server 8000 --bind 127.0.0.1
```

---

## Plotting Using PyGLVis

```bash
pip install glvis
```

To plot, please reference the `pyglivs.ipynb` file within the `plotting` folder. 
You will need to adjust the input files of `mesh` and the `GridFunction x` that you are plotting. 

---

## Core BESFEM Equations

![Alt text](equations.png)

**Mass Matrix:** used for any term with a time derivative.
```bash
AddDomainIntegrator(new mfem::MassIntegrator());
```

**Stiffness Matrix:** used for PDE terms involving ∇.
```bash
AddDomainIntegrator(new mfem::DiffusionIntegrator());
```

**Linear System Assembly:** Produces `A * X = B`
```bash
FormLinearSystem(ess_tdof_list, x, b, A, X, B);
```








