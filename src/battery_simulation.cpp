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

    cfg.init_cathode_particles = {0.35, 0.4, 0.30};
    cfg.init_anode_particles   = {0.2, 0.15, 0.10}; 

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

        const int np = static_cast<int>(state.anode_particles.size());

        if (mfem::Mpi::WorldRank() == 0)
        {
            for (int j = 0; j < np; ++j)
            {
                const double Xfr = state.anode_particles[j].concentration->GetLithiation();
                std::cout << "Xfr_" << j << " = " << Xfr << std::endl;
            }
        }

        // double VCell = 0.0;

        // ============================================================================
        // ===============================  TIME STEP LOOP  ===========================
        // ============================================================================

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
                    const int np = static_cast<int>(state.anode_particles.size());
                    std::vector<double> global_currents(np, 0.0);

                    UpdateAnodePairChemicalPotentials(state, geometry, domain_parameters);

                    *state.Rxn_gf = 0.0;
                    for (int j = 0; j < np; ++j)
                    {
                        *state.anode_particles[j].Rx_src = *state.anode_particles[j].Rxn_gf;
                        *state.Rxn_gf += *state.anode_particles[j].Rxn_gf;

                        std::vector<ConcentrationBase::PairCoupling> pair_terms;
                        Pairs(state, geometry, domain_parameters, j, pair_terms, np, t);

                        state.anode_particles[j].concentration->UpdateConcentration(*state.anode_particles[j].Rx_src, *state.anode_particles[j].Cn_gf,
                            *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);

                    }

                    state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf,
                        *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});

                    if (t > 0 && t % 50 == 0) {
                        state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                    }

                    for (int j = 0; j < np; ++j)
                    {
                        state.anode_particles[j].potential->AssembleSystem(*state.anode_particles[j].Cn_gf, *domain_parameters.ps[j], *state.anode_particles[j].ph_gf);
                    }
                    state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                    // *state.phA_gf = Constants::init_BvA;
                    // *state.phE_gf = Constants::init_BvE;

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    for (int j = 0; j < np; ++j)
                    {
                        state.anode_particles[j].reaction->TableExchangeCurrentDensity(*state.anode_particles[j].Cn_gf, *domain_parameters.AvEs[j]);
                        // while loop
                        while (globalerror_P > 1e-6 || globalerror_E > 1e-6) {
                            state.anode_particles[j].reaction->ButlerVolmer(*state.anode_particles[j].Rxn_gf, *state.anode_particles[j].Cn_gf,*state.CnE_gf,
                                *state.phA_gf, *state.phE_gf, *domain_parameters.AvEs[j]);
                            state.anode_particles[j].potential->UpdatePotential(*state.anode_particles[j].Rxn_gf, *state.anode_particles[j].ph_gf, *domain_parameters.ps[j], globalerror_P);
                            state.electrolyte_potential->UpdatePotential(*state.anode_particles[j].Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);
                        }
                        // while loop
                        state.anode_particles[j].reaction->TotalReactionCurrent(*state.anode_particles[j].Rxn_gf, global_currents[j]);
                    }

                    for (int j = 0; j < np; ++j)
                    {
                        double sgn = std::copysign(1.0, domain_parameters.gTrgPs[j] - global_currents[j]);
                        double dV  = Constants::dt * Constants::Vsr0 * sgn;

                        state.electrolyte_potential->BvE += dV;
                        *state.phE_gf += dV;
                    }



                    // ============================================================================
                    // ===============================  PRINT STATEMENTS  =========================
                    // ============================================================================


                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::ofstream outfile("anode_currents_mp.txt", std::ios::app);
                        outfile << "timestep: " << t;

                        for (int j = 0; j < np; ++j)
                        {
                            outfile << ", Current_" << j << " = " << global_currents[j] << ", Target_" << j << " = " << domain_parameters.gTrgPs[j];
                        }
                        outfile << std::endl;
                    }

                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::ofstream outfile("anode_concentrations_mp.txt", std::ios::app);
                        outfile << "timestep: " << t << " [ANODE HALF-CELL]";

                        for (int j = 0; j < np; ++j)
                        {
                            const double Xfr = state.anode_particles[j].concentration->GetLithiation();
                            outfile << ", Xfr_" << j << " = " << Xfr ;
                        }

                        outfile << std::endl;
                    }

                    
                }
                else
                {
                    const int np = static_cast<int>(state.cathode_particles.size());
                    std::vector<double> global_currents(np, 0.0);

                    UpdateCathodePairChemicalPotentials(state, geometry, domain_parameters);

                    *state.Rxn_gf = 0.0;
                    for (int j = 0; j < np; ++j)
                    {
                        *state.cathode_particles[j].Rx_src = *state.cathode_particles[j].Rxn_gf;
                        *state.Rxn_gf += *state.cathode_particles[j].Rxn_gf;
                        
                        std::vector<ConcentrationBase::PairCoupling> pair_terms;
                        Pairs(state, geometry, domain_parameters, j, pair_terms, np, t);
                        
                        state.cathode_particles[j].concentration->UpdateConcentration(*state.cathode_particles[j].Rx_src, *state.cathode_particles[j].Cn_gf,
                            *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
                    }

                    state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf,
                        *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});
                        
                    if (t > 0 && t % 50 == 0) {
                        state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                    }   
                    
                    for (int j = 0; j < np; ++j)
                    {
                        state.cathode_particles[j].potential->AssembleSystem(*state.cathode_particles[j].Cn_gf, *domain_parameters.ps[j], *state.cathode_particles[j].ph_gf);
                    }
                    state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    for (int j = 0; j < np; ++j)
                    {
                        state.cathode_particles[j].reaction->ExchangeCurrentDensity(*state.cathode_particles[j].Cn_gf, *domain_parameters.AvEs[j]);
                        // while loop
                        while (globalerror_P > 1e-6 || globalerror_E > 1e-6) {
                            state.cathode_particles[j].reaction->ButlerVolmer(*state.cathode_particles[j].Rxn_gf, *state.cathode_particles[j].Cn_gf,*state.CnE_gf,
                                *state.phC_gf, *state.phE_gf, *domain_parameters.AvEs[j]);
                            state.cathode_particles[j].potential->UpdatePotential(*state.cathode_particles[j].Rxn_gf, *state.cathode_particles[j].ph_gf, *domain_parameters.ps[j], globalerror_P);
                            state.electrolyte_potential->UpdatePotential(*state.cathode_particles[j].Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);
                        }
                        // while loop
                        state.cathode_particles[j].reaction->TotalReactionCurrent(*state.cathode_particles[j].Rxn_gf, global_currents[j]);
                    }

                    for (int j = 0; j < np; ++j)
                    {
                        double sgn = std::copysign(1.0, domain_parameters.gTrgPs[j] - global_currents[j]);
                        double dV  = Constants::dt * Constants::Vsr0 * sgn;

                        state.electrolyte_potential->BvE += dV;
                        *state.phE_gf += dV;
                    }


                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::ofstream outfile("cathode_currents_mp.txt", std::ios::app);
                        outfile << "timestep: " << t;

                        for (int j = 0; j < np; ++j)
                        {
                            outfile << ", Current_" << j << " = " << global_currents[j] << ", Target_" << j << " = " << domain_parameters.gTrgPs[j];
                        }
                        outfile << std::endl;
                    }

                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::ofstream outfile("cathode_concentrations_mp.txt", std::ios::app);
                        outfile << "timestep: " << t << " [CATHODE HALF-CELL]";

                        for (int j = 0; j < np; ++j)
                        {
                            const double Xfr = state.cathode_particles[j].concentration->GetLithiation();
                            outfile << ", Xfr_" << j << " = " << Xfr ;
                        }

                        outfile << std::endl;
                    }
                    
                }

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    std::vector<mfem::ParGridFunction*> anode_cn_fields;
                    anode_cn_fields.reserve(state.anode_particles.size());

                    for (auto &p : state.anode_particles)
                    {
                        anode_cn_fields.push_back(p.Cn_gf.get());
                    }

                    Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters,
                        anode_cn_fields, state.anode_out, 100);
                }
                else
                {
                    std::vector<mfem::ParGridFunction*> cathode_cn_fields;
                    cathode_cn_fields.reserve(state.cathode_particles.size());

                    for (auto &p : state.cathode_particles)
                    {
                        cathode_cn_fields.push_back(p.Cn_gf.get());
                    }

                    Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters,
                        cathode_cn_fields, state.cathode_out, 100);
                }
            }

        
        
        // else
        // {
        //     RunFullCellSimulation(state, geometry, domain_parameters, bc, adjust, outdir, cfg);
        // }
        


        }
    }
}