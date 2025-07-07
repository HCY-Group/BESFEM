# BESFEM


## How to Run the Code
1. Clone the repository
2. `cd BESFEM`
3. `make` (you may need to update the Makefile to point to your MFEM installation)
4. `cd bin`
5. Run a simulation: `mpirun -np 4 ./simulation -t d`

## Additional Command-Line Options


- `-m MeshFileName` This will allow you to change the .mesh file 
- `-d DistanceFileName` This will allow you to change the distance .txt file
- `-o Order` This will allow you to change the polynomial degree order
- `-t MeshType` This will allow you to change the type of mesh (r = rectangle, d = disk, v = voxel, c = circle)