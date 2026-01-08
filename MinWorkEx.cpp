#include "mfem.hpp"
#include "mpi.h"

int main(int argc, char *argv[]) {

    std::cout << "starting example" << std::endl;

    // Initialize MPI for parallel processing and HYPRE for solver setup
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);




	// Initialize some empty classes to make sure they work
	mfem::GridFunction test_gf;
    std::cout << "GridFunction initialized successfully" << std::endl;
	mfem::ParGridFunction test_pgf;
    std::cout << "ParGridFunction initialized successfully" << std::endl;
	mfem::Mesh test_mesh;
    std::cout << "Mesh initialized successfully" << std::endl;
	mfem::FiniteElementSpace test_fes;
    std::cout << "FiniteElementSpace initialized successfully" << std::endl;
	
	// try to initialize a distance solver
	//mfem::common::DistanceSolver test_ds;
	//mfem::DistanceSolver test_ds;
   // std::cout << "DistanceSolver initialized successfully" << std::endl;




    // Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    std::cout << "ending example" << std::endl;

    return 0;
}
