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

        // Initialize Concentrations & Potentials
        std::unique_ptr<mfem::ParGridFunction> CnA_gf, phA_gf, CnA_gf_psi;
        std::unique_ptr<CnA> anode_concentration;
        std::unique_ptr<PotA> anode_potential;

        std::unique_ptr<mfem::ParGridFunction> CnC_gf, phC_gf, CnC_gf_psi;
        std::unique_ptr<CnC> cathode_concentration;
        std::unique_ptr<PotC> cathode_potential;

        std::unique_ptr<mfem::ParGridFunction> CnE_gf, phE_gf, CnE_gf_psi;
        std::unique_ptr<CnE> electrolyte_concentration;
        std::unique_ptr<PotE> electrolyte_potential;

        std::unique_ptr<mfem::ParGridFunction> CnP_together;

        CnA_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        CnC_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        CnE_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        CnP_together = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        // always initialize electrolyte concentration & potential
        electrolyte_concentration = std::make_unique<CnE>(geometry, domain_parameters, bc, cfg.mode);
        CnE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        electrolyte_concentration->SetupField(*CnE_gf, Constants::init_CnE, *domain_parameters.pse);

        electrolyte_potential = std::make_unique<PotE>(geometry, domain_parameters, bc, cfg.mode);
        phE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        electrolyte_potential->SetupField(*phE_gf, Constants::init_BvE, *domain_parameters.pse);

        if (cfg.mode == sim::CellMode::HALF) // HALF-CELL
        {
            if(cfg.half_electrode == sim::Electrode::ANODE) {

                anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
                CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_concentration->SetupField(*CnA_gf, Constants::init_CnA, *domain_parameters.psi);

                anode_potential = std::make_unique<PotA>(geometry, domain_parameters, bc);
                phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_potential->SetupField(*phA_gf, Constants::init_BvA, *domain_parameters.psi);

            } else { // HALF-CATHODE

                cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
                CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_concentration->SetupField(*CnC_gf, Constants::init_CnC, *domain_parameters.psi);

                cathode_potential = std::make_unique<PotC>(geometry, domain_parameters, bc);
                phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_potential->SetupField(*phC_gf, Constants::init_BvC, *domain_parameters.psi);

            }
        } 
        else { // FULL-CELL

                anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
                CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_concentration->SetupField(*CnA_gf, Constants::init_CnA, *domain_parameters.psA);

                anode_potential = std::make_unique<PotA>(geometry, domain_parameters, bc);
                phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_potential->SetupField(*phA_gf, Constants::init_BvA, *domain_parameters.psA);

                cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
                CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_concentration->SetupField(*CnC_gf, Constants::init_CnC, *domain_parameters.psC);

                cathode_potential = std::make_unique<PotC>(geometry, domain_parameters, bc);
                phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_potential->SetupField(*phC_gf, Constants::init_BvC, *domain_parameters.psC);
        }

        // Initialize Reaction
        std::unique_ptr<mfem::ParGridFunction> Rxn_gf;
        std::unique_ptr<mfem::ParGridFunction> RxC_gf, RxA_gf;
        std::unique_ptr<Reaction> reaction;

        reaction = std::make_unique<Reaction>(geometry, domain_parameters);
        Rxn_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        reaction->Initialize(*Rxn_gf, Constants::init_Rxn);

        if (cfg.mode == sim::CellMode::FULL) { // FULL-CELL
            RxC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            RxA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            reaction->Initialize(*RxC_gf, Constants::init_RxC);
            reaction->Initialize(*RxA_gf, Constants::init_RxA);
        }

        // Set initial global current and cell voltage values  
        double global_current = 0.0;
        double global_current_A = 0.0;
        double global_current_C = 0.0;
        double VCell = 0.0;

        // Main Simulation Loop

        // ============================================================================
        // =================  HALF-CELL TIME STEPPING  ================================
        // ============================================================================

        if (cfg.mode == sim::CellMode::HALF) {

            int t = 0;

            if (cfg.half_electrode == sim::Electrode::ANODE) {
                VCell = anode_potential->BvA - electrolyte_potential->BvE;
            } else {
                VCell = cathode_potential->BvC - electrolyte_potential->BvE;
            }

            // for (int t = 0; t < cfg.num_timesteps; ++t) {

            while (VCell > 2.5){

                if (cfg.half_electrode == sim::Electrode::ANODE) {
                    
                    // ============================================================================
                    // =================  ANODE HALF-CELL TIME STEPPING  ==========================
                    // ============================================================================

                    anode_concentration->UpdateConcentration(*Rxn_gf, *CnA_gf, *domain_parameters.psi);
                    electrolyte_concentration->UpdateConcentration(*Rxn_gf, *CnE_gf, *domain_parameters.pse);

                    if (t > 0 && t % 500 == 0){
                        electrolyte_concentration->SaltConservation(*CnE_gf, *domain_parameters.pse);
                    }

                    anode_potential->AssembleSystem(*CnA_gf, *domain_parameters.psi, *phA_gf);
                    electrolyte_potential->AssembleSystem(*CnE_gf, *domain_parameters.pse, *phE_gf);

                    reaction->TableExchangeCurrentDensity(*CnA_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    while (globalerror_P > 1.0e-8 || globalerror_E > 1.0e-8) {
                        reaction->ButlerVolmer(*Rxn_gf, *CnA_gf, *CnE_gf, *phA_gf, *phE_gf);
                        anode_potential->UpdatePotential(*Rxn_gf, *phA_gf, *domain_parameters.psi, globalerror_P);
                        electrolyte_potential->UpdatePotential(*Rxn_gf, *phE_gf, *domain_parameters.pse, globalerror_E);
                        
                    }

                } else {
                
                    // ============================================================================
                    // ================  CATHODE HALF-CELL TIME STEPPING  =========================
                    // ============================================================================   


                    cathode_concentration->UpdateConcentration(*Rxn_gf, *CnC_gf, *domain_parameters.psi);
                    electrolyte_concentration->UpdateConcentration(*Rxn_gf, *CnE_gf, *domain_parameters.pse);

                    if (t > 0 && t % 500 == 0){
                        electrolyte_concentration->SaltConservation(*CnE_gf, *domain_parameters.pse);
                    }

                    cathode_potential->AssembleSystem(*CnC_gf, *domain_parameters.psi, *phC_gf);
                    electrolyte_potential->AssembleSystem(*CnE_gf, *domain_parameters.pse, *phE_gf);

                    reaction->ExchangeCurrentDensity(*CnC_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential
            
                    int iter = 0; 

                    while ((globalerror_P > 1.0e-8 || globalerror_E > 1.0e-8)) {
                        reaction->ButlerVolmer(*Rxn_gf, *CnC_gf, *CnE_gf, *phC_gf, *phE_gf);
                        cathode_potential->UpdatePotential(*Rxn_gf, *phC_gf, *domain_parameters.psi, globalerror_P);
                        electrolyte_potential->UpdatePotential(*Rxn_gf, *phE_gf, *domain_parameters.pse, globalerror_E);

                        if (iter == 50 && mfem::Mpi::WorldRank() == 0) {
                            std::cout << "WARNING: Nonlinear loop hit max iterations at timestep "
                                    << t << " | globalerror_P = " << globalerror_P
                                    << ", globalerror_E = " << globalerror_E << std::endl;                            
                        }

                        if (iter > 500) {
                            if (mfem::Mpi::WorldRank() == 0) {
                                std::cout << "ERROR: Nonlinear loop runaway at timestep " << t << std::endl;
                            }
                            break;
                        }

                        iter++;
                    }

                }

                reaction->TotalReactionCurrent(*Rxn_gf, global_current);

                double sgn = copysign(1.0, domain_parameters.gTrgI - global_current);
                double dV = Constants::dt * Constants::Vsr0 * sgn;
                electrolyte_potential->BvE += dV; // Adjust electrolyte potential based on target current
                *phE_gf += dV; // Update the grid function for electrolyte potential

                if (cfg.half_electrode == sim::Electrode::ANODE) {
                    VCell = anode_potential->BvA - electrolyte_potential->BvE;
                } else {
                    VCell = cathode_potential->BvC - electrolyte_potential->BvE;
                }

                if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0) {

                    std::ofstream outfile("half_cell_output.txt", std::ios::app);

                    const double Xfr = half_is_anode ? anode_concentration->GetLithiation() : cathode_concentration->GetLithiation();

                    outfile << "timestep: " << t << (half_is_anode ? " [ANODE HALF-CELL]" : " [CATHODE HALF-CELL]")
                    << ", Xfr = " << Xfr << ", VCell = " << VCell << ", BvE = " << electrolyte_potential->BvE
                    << (half_is_anode ? ", BvA = " : ", BvC = ") << (half_is_anode ? anode_potential->BvA : cathode_potential->BvC)
                    << ", current = " << global_current << ", target current = " << domain_parameters.gTrgI << std::endl;

                    outfile.close(); 
                }


                if (cfg.half_electrode == sim::Electrode::ANODE) {
                    Utils::SaveSimulationSnapshot(t, outdir, geometry, domain_parameters, *phA_gf, *phE_gf, 
                    *CnA_gf, *CnE_gf, *CnA_gf_psi, *CnE_gf_psi, 1000); 
                } else {
                    Utils::SaveSimulationSnapshot(t, outdir, geometry, domain_parameters, *phC_gf, *phE_gf, 
                    *CnC_gf, *CnE_gf, *CnC_gf_psi, *CnE_gf_psi, 1000); 
                }
 
                t += 1;

            } 
        }

        // ============================================================================
        // =================  FULL-CELL TIME STEPPING  ================================
        // ============================================================================

        if(cfg.mode == sim::CellMode::FULL) {
            
            double XfrA = anode_concentration->GetLithiation();
            double XfrC = cathode_concentration->GetLithiation();

            int t = 0;

            for (int t = 0; t < cfg.num_timesteps; ++t) {
            // while (XfrC < 0.85) {

            VCell = Constants::init_BvC - Constants::init_BvA;

            // while (VCell > 2.5) {

                anode_concentration->UpdateConcentration(*RxA_gf, *CnA_gf, *domain_parameters.psA);
                cathode_concentration->UpdateConcentration(*RxC_gf, *CnC_gf, *domain_parameters.psC);
                electrolyte_concentration->UpdateConcentration(*RxC_gf, *RxA_gf, *CnE_gf, *domain_parameters.pse); // with two inputs

                if (t > 0 && t % 500 == 0){
                    electrolyte_concentration->SaltConservation(*CnE_gf, *domain_parameters.pse);
                }

                cathode_potential->AssembleSystem(*CnC_gf, *domain_parameters.psC, *phC_gf);
                anode_potential->AssembleSystem(*CnA_gf, *domain_parameters.psA, *phA_gf);
                electrolyte_potential->AssembleSystem(*CnE_gf, *domain_parameters.pse, *phE_gf);

                reaction->ExchangeCurrentDensity(*CnC_gf, *CnA_gf); // with two inputs

                double globalerror_C = 1.0; // Error for cathode potential
                double globalerror_A = 1.0; // Error for anode potential
                double globalerror_E = 1.0; // Error for electrolyte potential

                double intlp = 0.0;

                while (globalerror_C > 1.0e-8 || globalerror_A > 1.0e-8 || globalerror_E > 1.0e-8) {
                
                    reaction->ButlerVolmer(*Rxn_gf, *RxC_gf, *RxA_gf, *CnC_gf, *CnA_gf, *CnE_gf, *phC_gf, *phA_gf, *phE_gf); // 9 inputs
                    cathode_potential->UpdatePotential(*RxC_gf, *phC_gf, *domain_parameters.psC, globalerror_C);
                    anode_potential->UpdatePotential(*RxA_gf, *phA_gf, *domain_parameters.psA, globalerror_A);
                    electrolyte_potential->UpdatePotential(*RxC_gf, *RxA_gf, *phE_gf, *domain_parameters.pse, globalerror_E);
                 }

                reaction->TotalReactionCurrent(*RxA_gf, global_current_A);
                reaction->TotalReactionCurrent(*RxC_gf, global_current_C);

                adjust.AdjustConstantCurrent(global_current_A, global_current_C, *anode_potential, *cathode_potential, *phA_gf, *phC_gf, VCell);

                XfrA = anode_concentration->GetLithiation();
                XfrC = cathode_concentration->GetLithiation();

                // Synchronize lithiation across ranks
                double global_XfrC;
                MPI_Allreduce(&XfrC, &global_XfrC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                global_XfrC /= mfem::Mpi::WorldSize();
                XfrC = global_XfrC; 

                if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0) {

                    std::ofstream outfile("full_cell_output.txt", std::ios::app);

                    outfile  << "timestep: " << t << " [FULL-CELL]" << ", XfrA = " << XfrA << ", XfrC = " << XfrC
                            << ", Anode current = " << global_current_A << ", Cathode current = " << global_current_C
                            << ", VCell = " << VCell << ", Target Current = " << domain_parameters.gTrgI << std::endl;

                    outfile.close(); 
                }
            
                Utils::SaveSimulationSnapshot(t, outdir, geometry, domain_parameters, *phA_gf, *phC_gf, *phE_gf, 
                    *CnA_gf, *CnC_gf, *CnE_gf, *CnA_gf_psi, *CnC_gf_psi, *CnE_gf_psi, *CnP_together, 100);
                 

                t += 1;

            } // end of FULL-CELL while loop

        } // end of FULL-CELL

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