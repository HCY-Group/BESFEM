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
                    const int np = static_cast<int>(state.cathode_particles.size());
                    UpdateCathodePairChemicalPotentials(state, geometry, domain_parameters);

                    // Assign Reaction Source
                    for (int j = 0; j < np; ++j)
                    {
                        *state.cathode_particles[j].Rx_src = Constants::init_RxC;
                    }

                    for (int j = 0; j < np; ++j) {
                        std::vector<ConcentrationBase::PairCoupling> pair_terms;

                        Pairs(state, geometry, domain_parameters, j, pair_terms, np);
                        
                        // for (int k = 0; k < np; ++k)
                        // {
                        //     if (j == k) { continue; }

                        //     const int a = std::min(j, k);
                        //     const int b = std::max(j, k);

                        //     ConcentrationBase::PairCoupling pair;
                        //     pair.sum_part = state.sum_pairs[a][b].get();
                        //     pair.weight   = domain_parameters.WeightPairs[a][b].get();
                        //     pair.grad_psi = domain_parameters.AvP_Pairs[a][b].get();

                        //     if (j < k)
                        //     {
                        //         pair.mu_self = state.mu_pair_a[a][b].get();
                        //         pair.mu_nbr  = state.mu_pair_b[a][b].get();
                        //     }
                        //     else
                        //     {
                        //         pair.mu_self = state.mu_pair_b[a][b].get();
                        //         pair.mu_nbr  = state.mu_pair_a[a][b].get();
                        //     }

                        //     pair_terms.push_back(pair);

                        //     if (mfem::Mpi::WorldRank() == 0 && t == 1)
                        //     {
                        //         std::cout << "[DEBUG] Pair (j,k) = (" << j << "," << k << ")"
                        //                 << " | stored as (a,b) = (" << a << "," << b << ")"; 
                        //         std::cout << std::endl;
                        //     }

                        // }

                        state.cathode_particles[j].concentration->UpdateConcentration(*state.cathode_particles[j].Rx_src, *state.cathode_particles[j].Cn_gf,
                            *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
                    }


                    
                }

                if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
                {
                    std::ofstream outfile("trial.txt", std::ios::app);
                    outfile << "timestep: " << t << " [CATHODE HALF-CELL]";
                    const int np = static_cast<int>(state.cathode_particles.size());

                    for (int j = 0; j < np; ++j)
                    {
                        const double Xfr = state.cathode_particles[j].concentration->GetLithiation();
                        outfile << ", Xfr_" << j << " = " << Xfr;
                    }

                    outfile << std::endl;
                    outfile.close();
                }

                std::vector<mfem::ParGridFunction*> cathode_cn_fields;
                cathode_cn_fields.reserve(state.cathode_particles.size());

                for (auto &p : state.cathode_particles)
                {
                    cathode_cn_fields.push_back(p.Cn_gf.get());
                }

                Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters,
                    cathode_cn_fields, state.cathode_out, 1000);
            }
        }

        
        
        // else
        // {
        //     RunFullCellSimulation(state, geometry, domain_parameters, bc, adjust, outdir, cfg);
        // }
        


    }
}