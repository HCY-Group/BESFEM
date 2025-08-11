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

// Sample runs:  mpirun -np 2 simulation
//               mpirun -np 1 simulation -m ../Mesh_3x90_T3.mesh
//               mpirun -np 6 simulation -m ../Mesh_40x30_3.mesh
//               mpirun -np 4 simulation -m ../Code_2D/II_1_bin.tif -d ../Code_2D/dsF_p.txt


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

    // Parse command-line options from MFEM
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&dsF_file, "-d", "--distance", "Distance file to use.");
    args.AddOption(&order, "-o", "--order", "Finite element polynomial degree.");
    args.AddOption(&mesh_type, "-t", "--type", "Mesh type (r for rectangle, c for circle, v for voxel).");
    args.ParseCheck();

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
    mfem::ParGridFunction CnCH_gf(geometry.parfespace.get());
    particle_concentration.Initialize(CnCH_gf, 0.02, *domain_parameters.psi);  // initial value: 2.02d-2

    // Electrolyte concentration initialization with a grid function and initial value
    CnE electrolyte_concentration(geometry, domain_parameters);
    mfem::ParGridFunction CnE_gf(geometry.parfespace.get());
    electrolyte_concentration.Initialize(CnE_gf, 0.001005, *domain_parameters.pse); 

    // Initialize Potential Classes

    // Particle potential initialization with a grid function and initial value
    PotP particle_potential(geometry, domain_parameters);
    mfem::ParGridFunction phP_gf(geometry.parfespace.get());
    // particle_potential.Initialize(phP_gf, 2.9395);
    particle_potential.Initialize(phP_gf, -0.1, *domain_parameters.psi);

    // Electrolyte potential initialization with a grid function and initial value
    PotE electrolyte_potential(geometry, domain_parameters);
    mfem::ParGridFunction phE_gf(geometry.parfespace.get());
    // electrolyte_potential.Initialize(phE_gf, -1.0);
    electrolyte_potential.Initialize(phE_gf, -0.4686, *domain_parameters.pse);

    
    // Initialize Reaction Class

    // Initialize the reaction with an empty reaction grid function and initial value
    Reaction reaction(geometry, domain_parameters);
    mfem::ParGridFunction Rxn_gf(geometry.parfespace.get());
    reaction.Initialize(Rxn_gf, 0.0);
    // reaction.Initialize(Rxn_gf, 1e-8);

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
    for (int t = 0; t < 30000 + 1; ++t) {

        particle_concentration.TimeStep(Rxn_gf, CnCH_gf, *domain_parameters.psi);
        electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, *domain_parameters.pse);

        // Only call SaltConservation every 500 timesteps
        if (t > 0 && t % 500 == 0){
                electrolyte_concentration.SaltConservation(CnE_gf, *domain_parameters.pse);
        }

        particle_potential.TimeStep(CnCH_gf, *domain_parameters.psi, phP_gf);
        electrolyte_potential.TimeStep(CnE_gf, *domain_parameters.pse, phE_gf, *electrolyte_concentration.CeVn);

        reaction.ExchangeCurrentDensity(CnCH_gf);

        double globalerror_P = 1.0; // Error for particle potential
        double globalerror_E = 1.0; // Error for electrolyte potential
 
        while (globalerror_P > 1.0e-9 || globalerror_E > 1.0e-9) {
            // Update reaction rates using the Butler-Volmer equation
            reaction.ButlerVolmer(Rxn_gf, CnCH_gf, CnE_gf, phP_gf, phE_gf);

            particle_potential.Advance(Rxn_gf, phP_gf, *domain_parameters.psi, globalerror_P);
            electrolyte_potential.Advance(Rxn_gf, phE_gf, *domain_parameters.pse, globalerror_E);
        }

        reaction.TotalReactionCurrent(Rxn_gf, global_current);

        double sgn = copysign(1.0, domain_parameters.gTrgI - global_current);
        double dV = Constants::dt * Constants::Vsr * sgn;
        particle_potential.BvP -= dV; // Adjust particle potential based on target current
        phP_gf -= dV; // Update the grid function for particle potential

        VCell = particle_potential.BvP - electrolyte_potential.BvE;

        // // t += 1;

        if (t % 200 == 0 && mfem::Mpi::WorldRank() == 0) {
            std::cout << "timestep: " << t
                    << ", Xfr = " << particle_concentration.GetLithiation()
                    << ", VCell = " << VCell << ", BvP = " << particle_potential.BvP
                    << std::endl;
        }

    }

    // Multiply Grid Functions for Error Calculations
    CnCH_gf *= *domain_parameters.psi;
    phP_gf *= *domain_parameters.psi;

    CnE_gf *= *domain_parameters.pse;
    phE_gf *= *domain_parameters.pse;

    // Save simulation outputs
    geometry.parallelMesh->Save("../outputs/Results/pmesh");
    domain_parameters.psi->Save("../outputs/Results/psi");
    domain_parameters.pse->Save("../outputs/Results/pse");

    // // // // // CnP_gf.Save("Results/CnP");
    CnCH_gf.Save("../outputs/Results/CnCH");
    CnE_gf.Save("../outputs/Results/CnE");
    phP_gf.Save("../outputs/Results/phP");
    phE_gf.Save("../outputs/Results/phE");
    Rxn_gf.Save("../outputs/Results/Rxn");

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


