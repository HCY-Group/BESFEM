#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "Constants.hpp"
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    Initialize_Geometry geometry;

    // Initialize and setup the handler
    geometry.InitializeMesh(Constants::mesh_file, MPI_COMM_WORLD, Constants::order);
    // mesh.SetupBoundaryConditions();

    // // get the finite element space
    // std::shared_ptr<mfem::ParFiniteElementSpace> fespace = geometry.GetParFiniteElementSpace();

    // Domain_Parameters domain_parameters(geometry);
    // domain_parameters.SetupDomainParameters(fespace);

    MPI_Finalize();
    return 0;

}
