#include "../mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Constants.hpp"
#include <mpi.h>
#include <iostream>

// Sample runs:  mpirun -np 4 simulation
//               mpirun -np 4 simulation -m ../Mesh_3x90_T3.mesh
//               mpirun -np 4 simulation -m ../Mesh_40x30_3.mesh
//               mpirun -np 4 simulation -m ../II_1_bin.tif
//               mpirun -np 4 ex0p -m ../data/square-disc.mesh -o 2 // MFEM Example


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Default values
    const char* mesh_file = Constants::mesh_file;
    int order = Constants::order;

    // Parse command-line options from MFEM
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&order, "-o", "--order", "Finite element polynomial degree.");
    args.ParseCheck();
    
    // Initialize Mesh & Geometry
    Initialize_Geometry geometry;
    geometry.InitializeMesh(mesh_file, MPI_COMM_WORLD, order);
    geometry.SetupBoundaryConditions();

    MPI_Finalize();
    return 0;

}

