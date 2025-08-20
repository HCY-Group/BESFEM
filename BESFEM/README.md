# BESFEM

## Things to Consider for Code Review
- Ideas on handling boundary conditions for different geometries (where should the user specify?)
- Concentrations Base currently inherits SolverSteps, should this just be usage rather than inheritance?
- How should users make changes? Through terminal command lines, editing Constant.cpp file, GUI?




## How to Run the Code
1. Clone the repository
3. `make` (you may need to update the Makefile to point to your MFEM installation)
4. `cd bin`
5. Run a simulation: `mpirun -np 1 simulation -m ../inputs/disk_Mesh_80x80x6.mesh -d ../inputs/disk_dsF_81x81x7.txt -t d -n 6`

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
