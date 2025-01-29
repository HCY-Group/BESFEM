#include "mfem.hpp"
#include "MeshBase.hpp"
#include "Constants.hpp"
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    MeshBase mesh;

    // Initialize and setup the handler
    mesh.InitializeMesh(Constants::mesh_file, MPI_COMM_WORLD, Constants::order);
    // mesh.SetupBoundaryConditions();

    // // Perform calculations
    // mesh.CalculatePhaseFields();

    // // Print results
    // mesh.PrintResults();

    MPI_Finalize();
    return 0;

}
