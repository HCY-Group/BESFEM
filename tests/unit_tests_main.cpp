#define CATCH_CONFIG_RUNNER

#include "mfem.hpp"
#include "mpi.h"
#include "../include/BESFEM_All.hpp"
#include "catch.hpp"

#include <chrono>
#include <iostream>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <vector>



int main( int argc, char* argv[]) {

     // Start measuring the program execution time
    using namespace std::chrono;
    auto program_start = high_resolution_clock::now();

    // Initialize MPI for parallel processing and HYPRE for solver setup
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int result = Catch::Session().run();

    
    // Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    return result;
}

