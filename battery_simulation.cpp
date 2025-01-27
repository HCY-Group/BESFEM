#include "mfem.hpp"
#include "Mesh_Test.hpp"
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    Mesh_Test mesh;

    // Initialize and setup the handler
    mesh.Initialize();
    // mesh.SetupBoundaryConditions();

    // // Perform calculations
    // mesh.CalculatePhaseFields();

    // // Print results
    // mesh.PrintResults();

    MPI_Finalize();
    return 0;

}
