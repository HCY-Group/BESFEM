#include "mfem.hpp"
#include "mpi.h"

#include "../inputs/Constants.hpp"
#include "../include/SimTypes.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/Concentrations_Base.hpp"
#include "../include/Potentials_Base.hpp"
#include "../include/CnC.hpp"
#include "../include/CnE.hpp"
#include "../include/CnA.hpp"
#include "../include/PotC.hpp"
#include "../include/PotA.hpp"
#include "../include/PotE.hpp"
#include "../include/Reaction.hpp"
// #include "../include/Current.hpp"

#include <chrono>
#include <iostream>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <ctime>

// ============================================================================

namespace fs = std::filesystem;

// Build an output directory
static std::string BuildRunOutdir(const char* mesh_file, int num_steps)
{
    // timestamp (local time)
    auto now   = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &now_c);
#else
    localtime_r(&now_c, &tm);
#endif
    std::ostringstream ts;
    ts << std::put_time(&tm, "%Y%m%d_%H%M%S");

    // mesh base name without directory or extension
    std::string mesh_name = fs::path(mesh_file).stem().string();

    std::ostringstream od;
    od << "../outputs/Results/"
       << ts.str()
       << "__nsteps=" << num_steps
       << "__mesh=" << mesh_name;

    return od.str();
}

// ============================================================================

struct SimulationConfig {
    sim::CellMode mode = sim::CellMode::HALF;
    sim::Electrode half_electrode = sim::Electrode::ANODE;
    const char* mesh_file = Constants::mesh_file;
    const char* dsF_file_A = Constants::dsF_file_A; // anode ψ distance
    const char* dsF_file_C = Constants::dsF_file_C; // cathode ψ distance
    const char* mesh_type = nullptr;
    int order = Constants::order;
    int num_timesteps = 1000; // Default number of timesteps
};


// ============================================================================
// ============================================================================
// ============================================================================


int main(int argc, char *argv[]) {

    // Start measuring the program execution time
    using namespace std::chrono;
    auto program_start = high_resolution_clock::now();

    // Initialize MPI for parallel processing and HYPRE for solver setup
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    SimulationConfig cfg;
    
    {
        const char* mode  = "half";    // "half" | "full"
        const char* half_elec  = "anode";   // "anode" | "cathode" (HALF only)
        const char* mesh_file  = cfg.mesh_file;
        const char* dsF_file_A   = Constants::dsF_file_A;
        const char* dsF_file_C   = Constants::dsF_file_C;
        const char* mesh_type  = nullptr;
        int         order      = cfg.order;
        int         num_timesteps     = cfg.num_timesteps;

        // Parse command-line options from MFEM
        mfem::OptionsParser args(argc, argv);
        args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
        args.AddOption(&dsF_file_C, "-dC", "--cathode-distance", "Cathode distance file to use.");
        args.AddOption(&dsF_file_A, "-dA", "--anode-distance", "Anode distance file to use.");
        args.AddOption(&order, "-o", "--order", "Finite element polynomial degree.");
        args.AddOption(&mesh_type, "-t", "--type", "Mesh type (r for rectangle, c for circle, v for voxel, d for disk).");
        args.AddOption(&num_timesteps, "-n", "--num-steps", "Number of timesteps to run the simulation.");
        args.AddOption(&mode,      "-mode", "--mode",      "Cell mode: half | full.");
        args.AddOption(&half_elec, "-elec", "--electrode", "HALF mode only: anode | cathode.");
        args.ParseCheck();

        cfg.mode        = (std::strcmp(mode, "full") == 0) ? sim::CellMode::FULL : sim::CellMode::HALF;
        cfg.half_electrode = (std::strcmp(half_elec, "cathode") == 0) ? sim::Electrode::CATHODE : sim::Electrode::ANODE;
        cfg.mesh_file   = mesh_file;
        cfg.dsF_file_A  = dsF_file_A;   // may be nullptr if not passed
        cfg.dsF_file_C  = dsF_file_C;   // may be nullptr if not passed

        // Sanity checks
        auto flag_present = [&](std::initializer_list<const char*> names)
        {
            for (int i = 1; i < argc; ++i)
                for (auto n : names)
                    if (std::strcmp(argv[i], n) == 0) return true;
            return false;
        };
        const bool used_dA = flag_present({"-dA","--anode-distance"});
        const bool used_dC = flag_present({"-dC","--cathode-distance"});

        if (cfg.mode == sim::CellMode::FULL) {
            if (!used_dA || !used_dC || !*cfg.dsF_file_A || !*cfg.dsF_file_C)
                mfem::mfem_error("FULL mode requires both -dA <file> and -dC <file>.");
        } else {
            const char* active = (cfg.half_electrode == sim::Electrode::CATHODE)
                               ? cfg.dsF_file_C : cfg.dsF_file_A;
            if (!active || !*active)
                mfem::mfem_error("HALF mode requires -dA (anode) or -dC (cathode) with a valid file.");
        }

        if (mfem::Mpi::WorldRank() == 0) {
            if (cfg.mode == sim::CellMode::HALF) { // ANODE
                if (cfg.half_electrode == sim::Electrode::ANODE) {
                    if (!used_dA) mfem::mfem_error("HALF-ANODE requires --dA <file> (do not pass --dC).");
                    if (used_dC)  mfem::mfem_error("HALF-ANODE: you passed --dC (cathode). Use --dA instead.");
                    if (!cfg.dsF_file_A || !std::strlen(cfg.dsF_file_A))
                        mfem::mfem_error("HALF-ANODE: empty --dA value.");
                } else { // CATHODE
                    if (!used_dC) mfem::mfem_error("HALF-CATHODE requires --dC <file> (do not pass --dA).");
                    if (used_dA)  mfem::mfem_error("HALF-CATHODE: you passed --dA (anode). Use --dC instead.");
                    if (!cfg.dsF_file_C || !std::strlen(cfg.dsF_file_C))
                        mfem::mfem_error("HALF-CATHODE: empty --dC value.");
                }
            } else { // FULL
                if (!used_dA || !used_dC)
                    mfem::mfem_error("FULL mode requires both --dA <file> and --dC <file>.");
                if (!cfg.dsF_file_A || !std::strlen(cfg.dsF_file_A))
                    mfem::mfem_error("FULL: empty --dA value.");
                if (!cfg.dsF_file_C || !std::strlen(cfg.dsF_file_C))
                    mfem::mfem_error("FULL: empty --dC value.");
            }
        }

        const char* active_dsF = (cfg.half_electrode == sim::Electrode::CATHODE) ? cfg.dsF_file_C : cfg.dsF_file_A;

        // Create timestamped output folder
        std::string outdir = BuildRunOutdir(mesh_file, num_timesteps);
        if (mfem::Mpi::WorldRank() == 0) {
            fs::create_directories(outdir);
            // Optional: write run metadata
            std::ofstream meta(outdir + "/run.txt");
            meta << "mesh_file=" << mesh_file << "\n"
                    << "dsF_file=" << dsF_file_A << "\n"
                    << "num_steps=" << num_timesteps << "\n"
                    << "order=" << order << "\n"
                    << "procs=" << mfem::Mpi::WorldSize() << "\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);

        bool half_mode     = (cfg.mode == sim::CellMode::HALF);
        bool half_is_anode = (cfg.half_electrode == sim::Electrode::ANODE);


        // ============================================================================
        // ===============================  START SIMULATION  =========================
        // ============================================================================

        // Initialize Mesh & Geometry
        Initialize_Geometry geometry;
        if (cfg.mode == sim::CellMode::HALF) {
            geometry.InitializeMesh(cfg.mesh_file, active_dsF, MPI_COMM_WORLD, cfg.order);
            geometry.SetupBoundaryConditions(sim::CellMode::HALF, cfg.half_electrode);
        } else {
            geometry.InitializeMesh(cfg.mesh_file, cfg.dsF_file_A, cfg.dsF_file_C, MPI_COMM_WORLD, cfg.order);
            geometry.SetupBoundaryConditions(sim::CellMode::FULL, sim::Electrode::BOTH);
        }

        geometry.parallelMesh->Save((outdir + "/pmesh").c_str());

        // Initialize and Calculate Domain Parameters
        Domain_Parameters domain_parameters(geometry);
        domain_parameters.SetupDomainParameters(mesh_type);

        // Initialize Concentrations & Potentials
        std::unique_ptr<mfem::ParGridFunction> CnA_gf, phA_gf;
        std::unique_ptr<CnA> anode_concentration;
        std::unique_ptr<PotA> anode_potential;

        std::unique_ptr<mfem::ParGridFunction> CnC_gf, phC_gf;
        std::unique_ptr<CnC> cathode_concentration;
        std::unique_ptr<PotC> cathode_potential;

        std::unique_ptr<mfem::ParGridFunction> CnE_gf, phE_gf;
        std::unique_ptr<CnE> electrolyte_concentration;
        std::unique_ptr<PotE> electrolyte_potential;


        // always initialize electrolyte concentration & potential
        electrolyte_concentration = std::make_unique<CnE>(geometry, domain_parameters);
        CnE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        electrolyte_concentration->Initialize(*CnE_gf, Constants::init_CnE, *domain_parameters.pse);

        electrolyte_potential = std::make_unique<PotE>(geometry, domain_parameters);
        phE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        electrolyte_potential->Initialize(*phE_gf, Constants::init_BvE, *domain_parameters.pse);

        if (cfg.mode == sim::CellMode::HALF) // HALF-CELL
        {
            if(cfg.half_electrode == sim::Electrode::ANODE) {

                anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
                CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_concentration->Initialize(*CnA_gf, Constants::init_CnA, *domain_parameters.psi);

                anode_potential = std::make_unique<PotA>(geometry, domain_parameters);
                phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_potential->Initialize(*phA_gf, Constants::init_BvA, *domain_parameters.psi);

            } else { // HALF-CATHODE

                cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
                CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_concentration->Initialize(*CnC_gf, Constants::init_CnC, *domain_parameters.psi);

                cathode_potential = std::make_unique<PotC>(geometry, domain_parameters);
                phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_potential->Initialize(*phC_gf, Constants::init_BvC, *domain_parameters.psi);

            }
        } 
        else { // FULL-CELL

                anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
                CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_concentration->Initialize(*CnA_gf, Constants::init_CnA, *domain_parameters.psA);

                anode_potential = std::make_unique<PotA>(geometry, domain_parameters);
                phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                anode_potential->Initialize(*phA_gf, Constants::init_BvA, *domain_parameters.psA);

                cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
                CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_concentration->Initialize(*CnC_gf, Constants::init_CnC, *domain_parameters.psC);

                cathode_potential = std::make_unique<PotC>(geometry, domain_parameters);
                phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                cathode_potential->Initialize(*phC_gf, Constants::init_BvC, *domain_parameters.psC);
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

        int t = 0;
        
        // Perform simulation over time steps
        // while (VCell > Constants::VCut) {
        // while (particle_concentration.GetLithiation() < 0.98) {

        // ============================================================================
        // ===============================  HALF-CELL  ================================
        // ============================================================================

        if (cfg.mode == sim::CellMode::HALF) {
            for (int t = 0; t < num_timesteps; ++t) {

                if (cfg.half_electrode == sim::Electrode::ANODE) {
                    
                    // ============================================================================
                    // ===============================  ANODE HALF-CELL  ==========================
                    // ============================================================================

                    anode_concentration->TimeStep(*Rxn_gf, *CnA_gf, *domain_parameters.psi);
                    electrolyte_concentration->TimeStep(*Rxn_gf, *CnE_gf, *domain_parameters.pse);

                    if (t > 0 && t % 100 == 0){
                        electrolyte_concentration->SaltConservation(*CnE_gf, *domain_parameters.pse);
                    }

                    anode_potential->TimeStep(*CnA_gf, *domain_parameters.psi, *phA_gf);
                    electrolyte_potential->TimeStep(*CnE_gf, *domain_parameters.pse, *phE_gf, electrolyte_concentration->CeVn);

                    reaction->TableExchangeCurrentDensity(*CnA_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    while (globalerror_P > 1.0e-8 || globalerror_E > 1.0e-8) {
                        reaction->ButlerVolmer(*Rxn_gf, *CnA_gf, *CnE_gf, *phA_gf, *phE_gf);
                        anode_potential->Advance(*Rxn_gf, *phA_gf, *domain_parameters.psi, globalerror_P);
                        electrolyte_potential->Advance(*Rxn_gf, *phE_gf, *domain_parameters.pse, globalerror_E);
                    }

                } else {
                
                    // ============================================================================
                    // ==============================  CATHODE HALF-CELL  =========================
                    // ============================================================================    

                    cathode_concentration->TimeStep(*Rxn_gf, *CnC_gf, *domain_parameters.psi);
                    electrolyte_concentration->TimeStep(*Rxn_gf, *CnE_gf, *domain_parameters.pse);

                    if (t > 0 && t % 100 == 0){
                        electrolyte_concentration->SaltConservation(*CnE_gf, *domain_parameters.pse);
                    }

                    cathode_potential->TimeStep(*CnC_gf, *domain_parameters.psi, *phC_gf);
                    electrolyte_potential->TimeStep(*CnE_gf, *domain_parameters.pse, *phE_gf, electrolyte_concentration->CeVn);

                    reaction->ExchangeCurrentDensity(*CnC_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential
            
                    while (globalerror_P > 1.0e-8 || globalerror_E > 1.0e-8) {
                        reaction->ButlerVolmer(*Rxn_gf, *CnC_gf, *CnE_gf, *phC_gf, *phE_gf);
                        cathode_potential->Advance(*Rxn_gf, *phC_gf, *domain_parameters.psi, globalerror_P);
                        electrolyte_potential->Advance(*Rxn_gf, *phE_gf, *domain_parameters.pse, globalerror_E);
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

                if (t % 5 == 0 && mfem::Mpi::WorldRank() == 0) {

                    const double Xfr = half_is_anode ? anode_concentration->GetLithiation()
                                                : cathode_concentration->GetLithiation();

                    std::cout << "timestep: " << t
                    << (half_is_anode ? " [ANODE HALF-CELL]" : " [CATHODE HALF-CELL]")
                    << ", Xfr = " << Xfr
                    << ", VCell = " << VCell
                    << ", BvE = " << electrolyte_potential->BvE
                    << (half_is_anode ? ", BvA = " : ", BvC = ")
                    << (half_is_anode ? anode_potential->BvA : cathode_potential->BvC)
                    << ", current = " << global_current
                    << std::endl;
                }
            }
        }

        // ============================================================================
        // ===============================  FULL-CELL  ================================
        // ============================================================================

        if(cfg.mode == sim::CellMode::FULL) {
            for (int t = 0; t < num_timesteps; ++t) {

                anode_concentration->TimeStep(*RxA_gf, *CnA_gf, *domain_parameters.psA);
                cathode_concentration->TimeStep(*RxC_gf, *CnC_gf, *domain_parameters.psC);
                electrolyte_concentration->TimeStep(*RxC_gf, *RxA_gf, *CnE_gf, *domain_parameters.pse); // with two inputs

                if (t > 0 && t % 500 == 0){
                    electrolyte_concentration->SaltConservation(*CnE_gf, *domain_parameters.pse);
                }

                cathode_potential->TimeStep(*CnC_gf, *domain_parameters.psC, *phC_gf);
                anode_potential->TimeStep(*CnA_gf, *domain_parameters.psA, *phA_gf);
                electrolyte_potential->TimeStep(*CnE_gf, *domain_parameters.pse, *phE_gf, electrolyte_concentration->CeVn);

                // electrolyte_potential->TimeStep(*CnE_gf, *domain_parameters.pse, *phE_gf, electrolyte_concentration->CeVn, cathode_concentration->Fct);

                reaction->ExchangeCurrentDensity(*CnC_gf, *CnA_gf); // with two inputs

                double globalerror_C = 1.0; // Error for cathode potential
                double globalerror_A = 1.0; // Error for anode potential
                double globalerror_E = 1.0; // Error for electrolyte potential

                double intlp = 0.0;

                // while (globalerror_C > 1.0e-8 || globalerror_A > 1.0e-8 || globalerror_E > 1.0e-8) {
                
                    reaction->ButlerVolmer(*Rxn_gf, *RxC_gf, *RxA_gf, *CnC_gf, *CnA_gf, *CnE_gf, *phC_gf, *phA_gf, *phE_gf); // 9 inputs
                    cathode_potential->Advance(*RxC_gf, *phC_gf, *domain_parameters.psC, globalerror_C);
                    anode_potential->Advance(*RxA_gf, *phA_gf, *domain_parameters.psA, globalerror_A);
                    electrolyte_potential->Advance(*RxC_gf, *RxA_gf, *phE_gf, *domain_parameters.pse, globalerror_E);
 
                    // std::cout << "  iter err: " << globalerror_C << ", " << globalerror_A << ", " << globalerror_E << std::endl;
                // }

                reaction->TotalReactionCurrent(*RxA_gf, global_current_A);
                reaction->TotalReactionCurrent(*RxC_gf, global_current_C);

                double Vsr;
                double dCrnt = abs(global_current_A - domain_parameters.gTrgI);
                if (dCrnt < abs(domain_parameters.gTrgI)*0.05) {Vsr = 0.025 * Constants::Vsr0;}
                else if (dCrnt < abs(domain_parameters.gTrgI)*0.10) {Vsr = 0.25 * Constants::Vsr0;}
                else {Vsr = 1.0 * Constants::Vsr0;}

                double sgnA = copysign(1.0, domain_parameters.gTrgI - abs(global_current_A));
                double dV_A = Constants::dt * Vsr * sgnA * 0.5;
                anode_potential->BvA += dV_A; // Adjust anode potential based on target current
                *phA_gf += dV_A; // Update the grid function for anode potential

                double sgnC = copysign(1.0, domain_parameters.gTrgI - global_current_C);
                double dV_C = Constants::dt * Vsr * sgnC * 2.0;
                cathode_potential->BvC -= dV_C; // Adjust cathode potential based on target current
                *phC_gf -= dV_C; // Update the grid function for cathode potential

                VCell = anode_potential->BvA - cathode_potential->BvC;


                if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0) {

                    const double XfrA = anode_concentration->GetLithiation();
                    const double XfrC = cathode_concentration->GetLithiation();

                    std::cout << "timestep: " << t << (" [FULL-CELL]") << ", XfrA = " << XfrA << ", XfrC = " << XfrC
                    << ", Anode current = " << global_current_A << ", Cathode current = " << global_current_C << ", VCell = " << VCell << ", Target Current = " << domain_parameters.gTrgI
                    << std::endl;
                }

            }
        }



        // ============================================================================
        // ===============================  SAVE OUTPUTS  =============================
        // ============================================================================

        // ============================================================================
        if (cfg.mode == sim::CellMode::HALF) {
            if (cfg.half_electrode == sim::Electrode::ANODE) {
                phA_gf->Save((outdir + "/phA_final").c_str());
                *CnA_gf *= *domain_parameters.psi;
                CnA_gf->Save((outdir + "/CnA_final").c_str());
            } else {
                phC_gf->Save((outdir + "/phC_final").c_str());
                *CnC_gf *= *domain_parameters.psi;
                CnC_gf->Save((outdir + "/CnC_final").c_str());
            }
        } else { // FULL
            phA_gf->Save((outdir + "/phA_final").c_str());
            phC_gf->Save((outdir + "/phC_final").c_str());
            phE_gf->Save((outdir + "/phE_final").c_str());
            CnA_gf->Save((outdir + "/CnA_final_raw").c_str());
            CnC_gf->Save((outdir + "/CnC_final_raw").c_str());
            CnE_gf->Save((outdir + "/CnE_final_raw").c_str());
            *CnA_gf *= *domain_parameters.psA;  CnA_gf->Save((outdir + "/CnA_final").c_str());
            *CnC_gf *= *domain_parameters.psC;  CnC_gf->Save((outdir + "/CnC_final").c_str());
            *CnE_gf *= *domain_parameters.pse;  CnE_gf->Save((outdir + "/CnE_final").c_str());
            *CnA_gf += *CnC_gf;  CnA_gf->Save((outdir + "/CnP_final").c_str()); // CnP = CnA + CnC
            RxA_gf->Save((outdir + "/RxA_final_raw").c_str());
            RxC_gf->Save((outdir + "/RxC_final_raw").c_str());
        }
    }

    std::cout << "Simulation complete." << std::endl;

    // Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    // End timing and output the total program execution time
    auto program_end = high_resolution_clock::now();
    std::cout << "Total Program Time: " 
              << duration_cast<seconds>(program_end - program_start).count() 
              << " seconds" << std::endl;

    return 0;
}