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
#include <vector>

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

    cfg.init_cathode_particles = {0.15, 0.20, 0.10};
    cfg.init_anode_particles   = {0.8, 0.7, 0.9}; 

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
            geometry.InitializeMesh(cfg.mesh_file, active_dsF, cfg.mesh_type, MPI_COMM_WORLD, cfg.order);
        } else {
            geometry.InitializeMesh(cfg.mesh_file, cfg.dsF_file_A, cfg.dsF_file_C, cfg.mesh_type, MPI_COMM_WORLD, cfg.order);
        }

        // Initialize and Calculate Domain Parameters
        Domain_Parameters domain_parameters(geometry);
        domain_parameters.SetupDomainParameters(cfg.mesh_type);

        // Initialize Boundary Conditions 
        BoundaryConditions bc(geometry, domain_parameters);
        if (cfg.mode == sim::CellMode::HALF) {
            bc.SetupBoundaryConditions(sim::CellMode::HALF, cfg.half_electrode);
        } else {
            bc.SetupBoundaryConditions(sim::CellMode::FULL, sim::Electrode::BOTH);
        }

        // Define Adjuster for Surface Voltage & Current
        Adjust adjust(geometry, domain_parameters);

        // Initialize Concentration & Potential & Reaction Fields
        SimulationState state;
        InitializeFields(state, geometry, domain_parameters, bc, cfg);

        // double VCell = 0.0;

        // Main time-stepping loop

        if (cfg.mode == sim::CellMode::HALF)
        {
            
            int t = 0;

            // if (cfg.half_electrode == sim::Electrode::ANODE) {
            //     VCell = anode_potential->BvA - electrolyte_potential->BvE;
            // } else {
            //     VCell = cathode_potential->BvC - electrolyte_potential->BvE;
            // }

            for (int t = 0; t < cfg.num_timesteps; ++t) {

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    // RunHalfCellSimulation(state, geometry, domain_parameters, bc, adjust, outdir, cfg);
                }
                else
                {
                    // RunHalfCellSimulation(state, geometry, domain_parameters, bc, adjust, outdir, cfg);
                }
            }
        }
        
        // else
        // {
        //     RunFullCellSimulation(state, geometry, domain_parameters, bc, adjust, outdir, cfg);
        // }
        


    }
}