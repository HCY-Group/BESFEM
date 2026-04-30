#include "mfem.hpp"
#include "mpi.h"
#include "../include/BESFEM_All.hpp"

#include <chrono>
#include <iostream>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <ctime>

int main(int argc, char *argv[]) {

    // Start measuring the program execution time
    using namespace std::chrono;
    auto program_start = high_resolution_clock::now();

    // Initialize MPI for parallel processing and HYPRE for solver setup
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    {

    SimulationConfig cfg = ParseSimulationArgs(argc, argv);
    ValidateConfig(cfg, argc, argv);

    std::string outdir = Utils::BuildRunOutdir(cfg.mesh_file, cfg.num_timesteps);
    if (mfem::Mpi::WorldRank() == 0)
        std::filesystem::create_directories(outdir);
    
    const char* active_dsF = (cfg.half_electrode == sim::Electrode::CATHODE) ? cfg.dsF_file_C : cfg.dsF_file_A;

    MPI_Barrier(MPI_COMM_WORLD);


    bool half_mode     = (cfg.mode == sim::CellMode::HALF);
    bool half_is_anode = (cfg.half_electrode == sim::Electrode::ANODE);


    // ============================================================================
    // ===============================  START SIMULATION  =========================
    // ============================================================================

    // Initialize Mesh & Geometry
    Initialize_Geometry geometry;
    if (cfg.mode == sim::CellMode::HALF) {
        geometry.InitializeMesh(cfg.mesh_file, active_dsF, cfg.mesh_type, MPI_COMM_WORLD, cfg.order, cfg.half_electrode);
    } else {
        geometry.InitializeMesh(cfg.mesh_file, cfg.dsF_file_A, cfg.dsF_file_C, cfg.mesh_type, MPI_COMM_WORLD, cfg.order);
    }

    // Initialize and Calculate Domain Parameters
    Domain_Parameters domain_parameters(geometry);
    domain_parameters.SetupDomainParameters(cfg.mesh_type);
    
    }

    if (mfem::Mpi::WorldRank() == 0) { std::cout << "Simulation complete.\n"; }

    // End timing and output the total program execution time
    auto program_end = high_resolution_clock::now();
    if (mfem::Mpi::WorldRank() == 0) {std::cout << "Total Program Time: " 
            << duration_cast<seconds>(program_end - program_start).count() 
            << " seconds" << std::endl;}

    // Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    return 0;
}