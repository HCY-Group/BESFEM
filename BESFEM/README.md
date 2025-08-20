# BESFEM

## Things to Consider for Code Review
- Other ways to organize the code into classes?
- Somewhere I am doing unnecessary calculations because OOP is MUCH slower than the non-OOP version.
- Memory management with pointers: unique vs shared?
- Ideas on handling boundary conditions for different geometries (where should the user specify?)
- Concentrations Base currently inherits SolverSteps, should this just be usage rather than inheritance?
- How should users make changes? Through terminal command lines, editing Constant.cpp file, GUI?
- Best way to switch particle concentration between diffusion and Cahn-Hilliard?




## How to Run the Code
1. Clone the repository
3. `make` (you may need to update the Makefile to point to your MFEM installation)
4. `cd bin`
5. Run a disk Cahn-Hilliard simulation: `mpirun -np 1 simulation -m ../inputs/disk_Mesh_80x80x6.mesh -d ../inputs/disk_dsF_81x81x7.txt -t d -n 6`
6. Run a rectangle diffusion simulation: 
    - ⚠️ make sure you go to `Constants.cpp` file and uncomment the section for rectangle
    - ⚠️ change `particle_concentration` in `simulation.cpp` to start with `CnP`and NOT `CnCH`
    - ⚠️ change all of the `Initialize` lines in `simulation.cpp` to the line underneath
    - `mpirun -np 1 simulation -m ../inputs/Mesh_3x90_r.mesh -d ../inputs/dsF_3x90_r.txt -t r -n 3`



## Additional Command-Line Options


- `-m MeshFileName` This will allow you to change the .mesh file 
- `-d DistanceFileName` This will allow you to change the distance .txt file
- `-o Order` This will allow you to change the polynomial degree order
- `-t MeshType` This will allow you to change the type of mesh (r = rectangle, d = disk, v = voxel, c = circle)
- `-n NumberOfTimeSteps` This will allow you to change the number of timesteps that the simulation runs


## How To Get Doxygen

1. Run `doxygen Doxyfile`
2. cd `html`
3. To view online when using VSCode: `python3 -m http.server 8000 --bind 127.0.0.1`


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

