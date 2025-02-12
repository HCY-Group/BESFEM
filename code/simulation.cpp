#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Constants.hpp"
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    // Initialize Mesh & Geometry
    Initialize_Geometry geometry;
    geometry.InitializeMesh(Constants::mesh_file, MPI_COMM_WORLD, Constants::order);
    geometry.SetupBoundaryConditions();

    MPI_Finalize();
    return 0;

}
