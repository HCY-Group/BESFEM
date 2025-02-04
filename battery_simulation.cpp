#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "Constants.hpp"
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    // Initialize Mesh & Geometry
    Initialize_Geometry geometry;
    geometry.InitializeMesh(Constants::mesh_file, MPI_COMM_WORLD, Constants::order);
    // mesh.SetupBoundaryConditions();

    // Initialize and Calculate Domain Parameters (psi, pse, AvB, AvP)
    Domain_Parameters domain_parameters(geometry);
    domain_parameters.SetupDomainParameters();

    MPI_Finalize();
    return 0;

}
