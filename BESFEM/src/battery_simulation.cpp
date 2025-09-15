#include "mfem.hpp"
#include "mpi.h"

#include "../inputs/Constants.hpp"
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

// Sample run: uses disk mesh and distance function and runs for 250 timesteps
// mpirun -np 6 battery_simulation -m ../inputs/disk_Mesh_80x80x6.mesh -da ../inputs/disk_dsF_81x81x7.txt -t d -n 250

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


// ---- Save all outputs (raw + projected) in one place ------------------------
static void SaveSimulationOutputs(
    const std::string& outdir,
    Initialize_Geometry& geometry,
    Domain_Parameters& dp,
    const mfem::ParGridFunction& CnP_gf,
    const mfem::ParGridFunction& CnE_gf,
    const mfem::ParGridFunction& phP_gf,
    const mfem::ParGridFunction& phE_gf,
    const mfem::ParGridFunction& Rxn_gf)
{
    // Ensure the directory exists before any rank writes into it
    if (mfem::Mpi::WorldRank() == 0) {
        std::filesystem::create_directories(outdir);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // 1) Save raw fields
    geometry.parallelMesh->Save((outdir + "/pmesh").c_str());
    dp.psi->Save((outdir + "/psi").c_str());
    dp.pse->Save((outdir + "/pse").c_str());
    dp.AvB->Save((outdir + "/AvB").c_str());
    dp.AvP->Save((outdir + "/AvP").c_str());

    CnP_gf.Save((outdir + "/CnP").c_str());
    CnE_gf.Save((outdir + "/CnE").c_str());
    phP_gf.Save((outdir + "/phP").c_str());
    phE_gf.Save((outdir + "/phE").c_str());
    Rxn_gf.Save((outdir + "/Rxn").c_str());

    // 2) CnP, phP are in the particle domain (mask by psi)
    {
        mfem::ParGridFunction pCnP(CnP_gf.ParFESpace()); pCnP = CnP_gf;
        pCnP *= *dp.psi;  pCnP.Save((outdir + "/pCnP").c_str());

        mfem::ParGridFunction pphP(phP_gf.ParFESpace()); pphP = phP_gf;
        pphP *= *dp.psi;  pphP.Save((outdir + "/pphP").c_str());
    }
    // CnE, phE are in the electrolyte domain (mask by pse)
    {
        mfem::ParGridFunction pCnE(CnE_gf.ParFESpace()); pCnE = CnE_gf;
        pCnE *= *dp.pse;  pCnE.Save((outdir + "/pCnE").c_str());

        mfem::ParGridFunction pphE(phE_gf.ParFESpace()); pphE = phE_gf;
        pphE *= *dp.pse;  pphE.Save((outdir + "/pphE").c_str());
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

enum class CellMode   { HALF, FULL };
enum class Electrode  { ANODE, CATHODE, BOTH }; // BOTH used only for FULL

struct SimulationConfig {
    CellMode    mode = CellMode::HALF;
    Electrode   half_electrode = Electrode::ANODE; // which one in HALF
    const char* mesh_file = Constants::mesh_file;
    const char* dsF_file_A = Constants::dsF_file_A; // anode ψ distance
    const char* dsF_file_C = Constants::dsF_file_C; // cathode ψ distance
    const char* mesh_type = nullptr;
    int order = Constants::order;
    int num_timesteps = 1000; // Default number of timesteps
};



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
        args.AddOption(&dsF_file_C, "-dC", "--cathode distance", "Cathode distance file to use.");
        args.AddOption(&dsF_file_A, "-dA", "--anode distance", "Anode distance file to use.");
        args.AddOption(&order, "-o", "--order", "Finite element polynomial degree.");
        args.AddOption(&mesh_type, "-t", "--type", "Mesh type (r for rectangle, c for circle, v for voxel, d for disk).");
        args.AddOption(&num_timesteps, "-n", "--num-steps", "Number of timesteps to run the simulation.");
        args.AddOption(&mode,      "-mode", "--mode",      "Cell mode: half | full.");
        args.AddOption(&half_elec, "-elec", "--electrode", "HALF mode only: anode | cathode.");
        args.ParseCheck();

        cfg.mode        = (std::strcmp(mode, "full") == 0) ? CellMode::FULL : CellMode::HALF;
        cfg.half_electrode = (std::strcmp(half_elec, "cathode") == 0) ? Electrode::CATHODE : Electrode::ANODE;
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

        const bool used_dA = flag_present({"-dA","--anode distance"});
        const bool used_dC = flag_present({"-dC","--cathode distance"});

        const char* active_dsF = nullptr;
        if (cfg.mode == CellMode::HALF) {
            active_dsF = (cfg.half_electrode == Electrode::CATHODE) ? cfg.dsF_file_C : cfg.dsF_file_A;
            if (!(active_dsF && std::strlen(active_dsF))) {
                mfem::mfem_error("HALF mode: missing/empty distance file (use -da/--dsA or -dc/--dsC).");
            }
        }

        if (mfem::Mpi::WorldRank() == 0) {
            if (cfg.mode == CellMode::HALF) { // ANODE
                if (cfg.half_electrode == Electrode::ANODE) {
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

        // Initialize Mesh & Geometry
        Initialize_Geometry geometry;
        // geometry.InitializeMesh(mesh_file, dsF_file_C, MPI_COMM_WORLD, order);
        geometry.InitializeMesh(mesh_file, active_dsF, MPI_COMM_WORLD, cfg.order);
        geometry.SetupBoundaryConditions();

        // Initialize and Calculate Domain Parameters (psi, pse, AvB, AvP)
        Domain_Parameters domain_parameters(geometry);
        domain_parameters.SetupDomainParameters(mesh_type);
        
        // Cathode concentration initialization with a grid function and initial value
        CnC cathode_concentration(geometry, domain_parameters);
        mfem::ParGridFunction CnC_gf(geometry.parfespace.get());
        cathode_concentration.Initialize(CnC_gf, Constants::init_CnC, *domain_parameters.psi);

        // Anode concentration initialization with a grid function and initial value
        CnA anode_concentration(geometry, domain_parameters);
        mfem::ParGridFunction CnA_gf(geometry.parfespace.get());
        anode_concentration.Initialize(CnA_gf, Constants::init_CnA, *domain_parameters.psi);  // initial value: 2.02d-2

        // Electrolyte concentration initialization with a grid function and initial value
        CnE electrolyte_concentration(geometry, domain_parameters);
        mfem::ParGridFunction CnE_gf(geometry.parfespace.get());
        electrolyte_concentration.Initialize(CnE_gf, Constants::init_CnE, *domain_parameters.pse); 

        // Cathode potential initialization with a grid function and initial value
        PotC cathode_potential(geometry, domain_parameters);
        mfem::ParGridFunction phC_gf(geometry.parfespace.get());
        cathode_potential.Initialize(phC_gf, Constants::init_BvC, *domain_parameters.psi);

        // Anode potential initialization with a grid function and initial value
        PotA anode_potential(geometry, domain_parameters);
        mfem::ParGridFunction phA_gf(geometry.parfespace.get());
        anode_potential.Initialize(phA_gf, Constants::init_BvA, *domain_parameters.psi);

        // Electrolyte potential initialization with a grid function and initial value
        PotE electrolyte_potential(geometry, domain_parameters);
        mfem::ParGridFunction phE_gf(geometry.parfespace.get());
        electrolyte_potential.Initialize(phE_gf, Constants::init_BvE, *domain_parameters.pse);

        // Initialize the reaction with an empty reaction grid function and initial value
        Reaction reaction(geometry, domain_parameters);
        mfem::ParGridFunction Rxn_gf(geometry.parfespace.get());
        reaction.Initialize(Rxn_gf, 0.0);
        // reaction.Initialize(Rxn_gf, 1.0e-8);

        // Initialize Current Class

        // Create the Current class to control current based on particle potential
        // Current current(geometry, domain_parameters);

        bool half_mode     = (cfg.mode == CellMode::HALF);
        bool half_is_anode = (cfg.half_electrode == Electrode::ANODE);

        // Set initial global current and cell voltage values  
        double global_current = 0.0;
        double VCell = 0.0;

        // Main Simulation Loop

        int t = 0;
        
        // Perform simulation over time steps
        // while (VCell > Constants::VCut) {
        // while (particle_concentration.GetLithiation() < 0.98) {

        // ============================================================================
        // ===============================  HALF-CELL  ================================
        // ============================================================================

        if (half_mode) {
            for (int t = 0; t < num_timesteps; ++t) {

                if (half_is_anode) {
                    
                    // ============================================================================
                    // ===============================  ANODE HALF-CELL  ==========================
                    // ============================================================================

                    anode_concentration.TimeStep(Rxn_gf, CnA_gf, *domain_parameters.psi);
                    electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, *domain_parameters.pse);

                    if (t > 0 && t % 100 == 0){
                        electrolyte_concentration.SaltConservation(CnE_gf, *domain_parameters.pse);
                    }

                    anode_potential.TimeStep(CnA_gf, *domain_parameters.psi, phA_gf);
                    electrolyte_potential.TimeStep(CnE_gf, *domain_parameters.pse, phE_gf, electrolyte_concentration.CeVn);

                    reaction.TableExchangeCurrentDensity(CnA_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    mfem::StopWatch sw1;
                    sw1.Start();
            
                    while (globalerror_P > 1.0e-8 || globalerror_E > 1.0e-8) {

                        // Update reaction rates using the Butler-Volmer equation
                        reaction.ButlerVolmer(Rxn_gf, CnA_gf, CnE_gf, phA_gf, phE_gf);

                        anode_potential.Advance(Rxn_gf, phA_gf, *domain_parameters.psi, globalerror_P);
                        electrolyte_potential.Advance(Rxn_gf, phE_gf, *domain_parameters.pse, globalerror_E);
                    }

                    sw1.Stop();
                    if(mfem::Mpi::WorldRank() == 0) {
                        std::cout << "Timestep " << t << " advance time: " << sw1.RealTime() << " seconds" << std::endl;
                    }

                } else {
                
                    // ============================================================================
                    // ==============================  CATHODE HALF-CELL  =========================
                    // ============================================================================    

                    cathode_concentration.TimeStep(Rxn_gf, CnC_gf, *domain_parameters.psi);
                    electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, *domain_parameters.pse);

                    if (t > 0 && t % 100 == 0){
                        electrolyte_concentration.SaltConservation(CnE_gf, *domain_parameters.pse);
                    }

                    cathode_potential.TimeStep(CnC_gf, *domain_parameters.psi, phC_gf);
                    electrolyte_potential.TimeStep(CnE_gf, *domain_parameters.pse, phE_gf, electrolyte_concentration.CeVn);

                    reaction.ExchangeCurrentDensity(CnC_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    mfem::StopWatch sw1;
                    sw1.Start();
            
                    while (globalerror_P > 1.0e-8 || globalerror_E > 1.0e-8) {

                        // Update reaction rates using the Butler-Volmer equation
                        reaction.ButlerVolmer(Rxn_gf, CnC_gf, CnE_gf, phC_gf, phE_gf);

                        cathode_potential.Advance(Rxn_gf, phC_gf, *domain_parameters.psi, globalerror_P);
                        electrolyte_potential.Advance(Rxn_gf, phE_gf, *domain_parameters.pse, globalerror_E);
                    }

                    sw1.Stop();
                    if(mfem::Mpi::WorldRank() == 0) {
                        std::cout << "Timestep " << t << " advance time: " << sw1.RealTime() << " seconds" << std::endl;
                    }


                }

                reaction.TotalReactionCurrent(Rxn_gf, global_current);

                double sgn = copysign(1.0, domain_parameters.gTrgI - global_current);
                double dV = Constants::dt * Constants::Vsr * sgn;
                electrolyte_potential.BvE += dV; // Adjust electrolyte potential based on target current
                phE_gf += dV; // Update the grid function for electrolyte potential

                if (half_is_anode) {
                    VCell = anode_potential.BvA - electrolyte_potential.BvE;
                } else {
                    VCell = cathode_potential.BvC - electrolyte_potential.BvE;
                }

                if (t % 1 == 0 && mfem::Mpi::WorldRank() == 0) {

                    const double Xfr = half_is_anode ? anode_concentration.GetLithiation()
                                                : cathode_concentration.GetLithiation();

                    std::cout << "timestep: " << t
                    << (half_is_anode ? " [ANODE HALF-CELL]" : " [CATHODE HALF-CELL]")
                    << ", Xfr = " << Xfr
                    << ", VCell = " << VCell
                    << ", BvE = " << electrolyte_potential.BvE
                    << (half_is_anode ? ", BvA = " : ", BvC = ")
                    << (half_is_anode ? anode_potential.BvA : cathode_potential.BvC)
                    << ", current = " << global_current
                    << std::endl;
                        }
            }
            // Save final outputs
            if (half_is_anode) {
                SaveSimulationOutputs(outdir, geometry, domain_parameters, CnA_gf, CnE_gf, phA_gf, phE_gf, Rxn_gf);
            } else { 
                SaveSimulationOutputs(outdir, geometry, domain_parameters, CnC_gf, CnE_gf, phC_gf, phE_gf, Rxn_gf);
            }
        }
    }
    // Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    // End timing and output the total program execution time
    auto program_end = high_resolution_clock::now();
    std::cout << "Total Program Time: " 
              << duration_cast<seconds>(program_end - program_start).count() 
              << " seconds" << std::endl;

    // std::cout << "t skip value: " << t_skip << std::endl;

    return 0;
}




