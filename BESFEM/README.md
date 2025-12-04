# BESFEM

## Things to Consider for Code Review
- Other ways to organize the code into classes?
- Thoughts on boundary conditions class?
- How should users make changes? Through terminal command lines, editing Constant.cpp file, GUI?
- Best way to switch particle concentration between diffusion and Cahn-Hilliard? Right now, the anode uses Cahn-Hilliard



## How to Run the Code
1. Clone the repository
3. `make` (you may need to update the Makefile to point to your MFEM installation)
4. `cd bin`
5. Example for a full cell code: `mpirun -np 2 battery_simulation -mode full -m ../inputs/mesh/Mesh_40x60x3_3D_disk_full.mesh -dA ../inputs/distance/dsF_A_40x60x3_3D_disk_full.txt -dC ../inputs/distance/dsF_C_40x60x3_3D_disk_full.txt -t ml`
6. Example for a half cell (cathode) code: `mpirun -np 1 battery_simulation -mode half -elec cathode -m ../inputs/mesh/Mesh_40x60_F00.mesh  -dC ../inputs/distance/dsFC_41x61_F00.txt -t ml`


## Additional Command-Line Options


- `-m MeshFileName` This will allow you to change the .mesh file 
- `-mode CellMode` This will allow you to change if you are simulating a half or full cell. 
- `-elec Electrode` If you are doing a half cell simulation, this will allow you to specify which electrode you are simulating. 
- `-dC CathodeDistanceFileName` This will allow you to change the distance .txt file for the cathode
- `-dA AnodeDistanceFileName` This will allow you to change the distance .txt file for the anode
- `-o Order` This will allow you to change the polynomial degree order
- `-t MeshType` This will allow you to change the type of mesh (ml = MATLAB, v = voxel)
- `-n NumberOfTimeSteps` This will allow you to change the number of timesteps that the simulation runs


## How To Get Doxygen

1. Run `module load Doxygen`
2. Run `doxygen Doxyfile`
3. cd `html`
4. To view online when using VSCode: `python3 -m http.server 8000 --bind 127.0.0.1`


## How to Plot in Python (currently only supports serial)

1. Create a folder in your directory called `outputs`
2. Within `outputs`, create another folder called `Results`
3. When you run the code, the results will spit out inside another folder in `Results`
4. Edit the path to which you would like to plot in `mfem_plot.ipynb`:
    - `mesh_path = "outputs/Results/20250820_092154__nsteps=3__mesh=Mesh_3x90_r/pmesh.000000"`
    - `field_path = "outputs/Results/20250820_092154__nsteps=3__mesh=Mesh_3x90_r/pCnCH.000000"`
5. `Run All` and scroll down to the graph


## How to Remove Unwanted Results Files

1. `cd outputs/Results`
2. `rm -rf 2025*/`



## BESFEM Equations

![Alt text](equations.png)

Additional Notes:
- A mass matrix is used for anything with a time derivative
    - `AddDomainIntegrator(new mfem::MassIntegrator)`

- A stiffness matrix is used for anything with ∇ (PDEs)
    - `AddDomainIntegrator(new mfem::DiffusionIntegrator)`
    
- `FormLinearSystem(Array, x, b, A, X, B)` results in `A(X) = B`

