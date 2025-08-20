#include "mfem.hpp"
#include "mpi.h"

#include "../inputs/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/Concentrations_Base.hpp"
#include "../include/Potentials_Base.hpp"
#include "../include/CnP.hpp"
#include "../include/CnE.hpp"
#include "../include/CnCH.hpp"
#include "../include/PotP.hpp"
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

// Sample run: uses disk mesh and distance function and runs for 2000 timesteps
// mpirun -np 6 simulation -m ../inputs/disk_Mesh_80x80x6.mesh -d ../inputs/disk_dsF_81x81x7.txt -t d -n 250

// ============================================================================

namespace fs = std::filesystem;

// Build an output directory like:
// ../outputs/Results/20250812_0915__nsteps=50001__mesh=disk_Mesh_80x80x6_01
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
    // Default values
    const char* mesh_file = Constants::mesh_file;
    const char* dsF_file = Constants::dsF_file;
    const char* mesh_type = nullptr;
    int order = Constants::order;
    int num_timesteps = 1000; // Default number of timesteps
    int model_choice = 0; // Default model choice (0 for diffusion, 1 for Cahn-Hilliard)

    // Parse command-line options from MFEM
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&dsF_file, "-d", "--distance", "Distance file to use.");
    args.AddOption(&order, "-o", "--order", "Finite element polynomial degree.");
    args.AddOption(&mesh_type, "-t", "--type", "Mesh type (r for rectangle, c for circle, v for voxel, d for disk).");
    args.AddOption(&num_timesteps, "-n", "--num-steps", "Number of timesteps to run the simulation.");
    // args.AddOption(&model_choice, "-pm", "--model-choice", "Model choice for particle (0 for diffusion, 1 for Cahn-Hilliard).");
    args.ParseCheck();

    std::string mesh_type_str = (mesh_type ? mesh_type : "");

    // Create timestamped output folder
    std::string outdir = BuildRunOutdir(mesh_file, num_timesteps);
    if (mfem::Mpi::WorldRank() == 0) {
        fs::create_directories(outdir);
        // Optional: write run metadata
        std::ofstream meta(outdir + "/run.txt");
        meta << "mesh_file=" << mesh_file << "\n"
                << "dsF_file=" << dsF_file << "\n"
                << "num_steps=" << num_timesteps << "\n"
                << "order=" << order << "\n"
                << "procs=" << mfem::Mpi::WorldSize() << "\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize Mesh & Geometry
    Initialize_Geometry geometry;
    geometry.InitializeMesh(mesh_file, dsF_file, MPI_COMM_WORLD, order);
    geometry.SetupBoundaryConditions();

    // Initialize and Calculate Domain Parameters (psi, pse, AvB, AvP)
    Domain_Parameters domain_parameters(geometry);
    domain_parameters.SetupDomainParameters(mesh_type);

    // Initialize Concentration Classes
    
    // // Particle concentration initialization with a grid function and initial value
    // CnP particle_concentration(geometry, domain_parameters);
    // mfem::ParGridFunction CnP_gf(geometry.parfespace.get());
    // particle_concentration.Initialize(CnP_gf, 0.3, *domain_parameters.psi);

    // Particle concentration initialization with a grid function and initial value
    CnCH particle_concentration(geometry, domain_parameters);
    mfem::ParGridFunction CnP_gf(geometry.parfespace.get());
    particle_concentration.Initialize(CnP_gf, 2.02e-2, *domain_parameters.psi);  // initial value: 2.02d-2
    // particle_concentration.Initialize(CnP_gf, 0.3, *domain_parameters.psi);  // initial value: 2.02d-2


    // Electrolyte concentration initialization with a grid function and initial value
    CnE electrolyte_concentration(geometry, domain_parameters);
    mfem::ParGridFunction CnE_gf(geometry.parfespace.get());
    electrolyte_concentration.Initialize(CnE_gf, 0.001005, *domain_parameters.pse); 
    // electrolyte_concentration.Initialize(CnE_gf, 0.001, *domain_parameters.pse); // initial value: 0.001005

    // Initialize Potential Classes

    // Particle potential initialization with a grid function and initial value
    PotP particle_potential(geometry, domain_parameters);
    mfem::ParGridFunction phP_gf(geometry.parfespace.get());
    particle_potential.Initialize(phP_gf, -0.1, *domain_parameters.psi);
    // particle_potential.Initialize(phP_gf, 2.9395, *domain_parameters.psi);


    // Electrolyte potential initialization with a grid function and initial value
    PotE electrolyte_potential(geometry, domain_parameters);
    mfem::ParGridFunction phE_gf(geometry.parfespace.get());
    electrolyte_potential.Initialize(phE_gf, -0.4686, *domain_parameters.pse);
    // electrolyte_potential.Initialize(phE_gf, -1.0, *domain_parameters.pse);


    // Initialize Reaction Class

    // Initialize the reaction with an empty reaction grid function and initial value
    Reaction reaction(geometry, domain_parameters);
    mfem::ParGridFunction Rxn_gf(geometry.parfespace.get());
    reaction.Initialize(Rxn_gf, 0.0);
    // reaction.Initialize(Rxn_gf, 1.0e-8);

    // Initialize Current Class

    // Create the Current class to control current based on particle potential
    // Current current(geometry, domain_parameters);

    // Set initial global current and cell voltage values  
    double global_current = 0.0;
    double VCell = particle_potential.BvP - electrolyte_potential.BvE;

    // Main Simulation Loop

    int t = 0;
    
    // Perform simulation over time steps
    // while (VCell > Constants::VCut) {
    // while (particle_concentration.GetLithiation() < 0.98) {
    for (int t = 0; t < num_timesteps; ++t) {

        particle_concentration.TimeStep(Rxn_gf, CnP_gf, *domain_parameters.psi);
        electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, *domain_parameters.pse);

        // Only call SaltConservation every 100 timesteps
        if (t > 0 && t % 100 == 0){
                electrolyte_concentration.SaltConservation(CnE_gf, *domain_parameters.pse);
        }

        particle_potential.TimeStep(CnP_gf, *domain_parameters.psi, phP_gf);
        electrolyte_potential.TimeStep(CnE_gf, *domain_parameters.pse, phE_gf, electrolyte_concentration.CeVn);  
        
        if (mesh_type_str == "d") {
            reaction.TableExchangeCurrentDensity(CnP_gf);
        } else {
            reaction.ExchangeCurrentDensity(CnP_gf);
        }
        

        double globalerror_P = 1.0; // Error for particle potential
        double globalerror_E = 1.0; // Error for electrolyte potential
 
        while (globalerror_P > 1.0e-8 || globalerror_E > 1.0e-8) {
            // Update reaction rates using the Butler-Volmer equation
            reaction.ButlerVolmer(Rxn_gf, CnP_gf, CnE_gf, phP_gf, phE_gf);

            particle_potential.Advance(Rxn_gf, phP_gf, *domain_parameters.psi, globalerror_P);
            electrolyte_potential.Advance(Rxn_gf, phE_gf, *domain_parameters.pse, globalerror_E);
        }

        reaction.TotalReactionCurrent(Rxn_gf, global_current);

        double sgn = copysign(1.0, domain_parameters.gTrgI - global_current);
        double dV = Constants::dt * Constants::Vsr * sgn;
        electrolyte_potential.BvE += dV; // Adjust electrolyte potential based on target current
        phE_gf += dV; // Update the grid function for electrolyte potential

        VCell = particle_potential.BvP - electrolyte_potential.BvE;

        if (t % 200 == 0 && mfem::Mpi::WorldRank() == 0) {
            std::cout << "timestep: " << t
                    << ", Xfr = " << particle_concentration.GetLithiation()
                    << ", VCell = " << VCell << ", BvE = " << electrolyte_potential.BvE
                    << ", current = " << global_current
                    << std::endl;
        }

    }

    // Save outputs into the timestamped folder
    geometry.parallelMesh->Save((outdir + "/pmesh").c_str());
    domain_parameters.psi->Save((outdir + "/psi").c_str());
    domain_parameters.pse->Save((outdir + "/pse").c_str());
    domain_parameters.AvB->Save((outdir + "/AvB").c_str());
    domain_parameters.AvP->Save((outdir + "/AvP").c_str());
    
    CnP_gf.Save((outdir + "/CnP").c_str());
    CnE_gf.Save((outdir + "/CnE").c_str());
    phP_gf.Save((outdir + "/phP").c_str());
    phE_gf.Save((outdir + "/phE").c_str());
    Rxn_gf.Save((outdir + "/Rxn").c_str());

    // Multiply Grid Functions for Error Calculations
    CnP_gf *= *domain_parameters.psi;
    phP_gf *= *domain_parameters.psi;
    CnE_gf *= *domain_parameters.pse;
    phE_gf *= *domain_parameters.pse;

    CnP_gf.Save((outdir + "/pCnP").c_str());
    CnE_gf.Save((outdir + "/pCnE").c_str());
    phP_gf.Save((outdir + "/pphP").c_str());
    phE_gf.Save((outdir + "/pphE").c_str());

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




