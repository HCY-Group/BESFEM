#include "mfem.hpp"
#include "mpi.h"

#include "../code/Constants.hpp"
#include "../code/Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "Concentrations_Base.hpp"
#include "Potentials_Base.hpp"
#include "CnP.hpp"
#include "CnE.hpp"
#include "CnCH.hpp"
#include "PotP.hpp"
#include "PotE.hpp"
#include "Reaction.hpp"
// #include "Current.hpp"

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
    electrolyte_concentration.Initialize(CnE_gf, 0.001, *domain_parameters.pse); 

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

    // // // Initialize Current Class

    // // Create the Current class to control current based on particle potential
    // Current current(geometry, domain_parameters);

    // Set initial global current and cell voltage values  
    double global_current = 0.0;
    double VCell = particle_potential.BvP - electrolyte_potential.BvE;


    // Main Simulation Loop

    // int t = 0;
    // int local_nE = geometry.nE;  // Local number of elements on this MPI process
    // int global_nE = 0;

    // MPI_Allreduce(&local_nE, &global_nE, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // int t_skip = std::max(1, static_cast<int>(std::ceil(global_nE / 20.0)));
    
    // Perform simulation over time steps
    for (int t = 0; t < 1000 + 1; ++t) {
        particle_concentration.TimeStep(Rxn_gf, CnCH_gf, *domain_parameters.psi);
        electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, *domain_parameters.pse);

        if (t % 50 == 0 && mfem::Mpi::WorldRank() == 0) {
            std::cout << "timestep: " << t
                    << ", Xfr = " << particle_concentration.GetLithiation()
                    << ", VCell = " << VCell
                    << std::endl;
        }


    }
    // // while ( VCell > Constants::VCut) {
    // while ( Xfr < 0.98 ) {

    //     // if (rank == 0) {
    //     //     std::cout << "Before timestep: Cn norm = " << CnCH_gf.Norml2() << std::endl;
    //     // }
    //     // Step 1: Update concentrations for both particle and electrolyte phases
    //     // particle_concentration.TimeStep(Rxn_gf, CnP_gf, *domain_parameters.psi);
        //     // check_nan("CnCH", CnCH_gf, t);

    //     // electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, *domain_parameters.pse);
    //     // // check_nan("CnE", CnE_gf, t);

    //     // // // // // if (t % t_skip == 0) {
    //     // // // //     // Step 2: Update potentials for both particle and electrolyte phases
    //     //     // particle_potential.TimeStep(CnCH_gf, *domain_parameters.psi, phP_gf);
    //     //     // electrolyte_potential.TimeStep(CnE_gf, *domain_parameters.pse, phE_gf);
    //     // // // // // // }

    //     // // // // Step 3: Compute rate constants and exchange current density at the interface     
    //     // reaction.ExchangeCurrentDensity(CnCH_gf); 
    //     // // reaction.WriteKfwToFile("Kfw_output.txt");


    //     // // // Step 4: Iteratively solve for global reaction rates using Butler-Volmer kinetics
    //     // double globalerror_P = 1.0; // Error for particle potential
    //     // double globalerror_E = 1.0; // Error for electrolyte potential

    //     // // // // if (t % t_skip == 0) {
    //     //     // while (globalerror_P > 1.0e-9 || globalerror_E > 1.0e-9) {

    //     //         // Update reaction rates using the Butler-Volmer equation
    //             // reaction.ButlerVolmer(Rxn_gf, CnCH_gf, CnE_gf, phP_gf, phE_gf);
    //     //         // check_nan("Rxn", Rxn_gf, t);


    //     // // //         // Calculate global errors for particle and electrolyte potentials
    //     //         // particle_potential.CalculateGlobalError(Rxn_gf, phP_gf, *domain_parameters.psi, globalerror_P);
    //     //         // electrolyte_potential.CalculateGlobalError(Rxn_gf, phE_gf, *domain_parameters.pse, globalerror_E);
            
    //     //     // }
    //     // // // }

    //     // // // Step 5: Compute the total reaction current
    //     // reaction.TotalReactionCurrent(Rxn_gf, global_current);

    //     // // // Step 6: Adjust the particle potential based on the target constant current
    //     // current.Constant(phP_gf, global_current);

    //     // Step 7: Update the cell voltage based on the particle and electrolyte potentials
    //     VCell = BvP - BvE;
        
    //     if (t % 200 == 0) {
    //         std::cout << "timestep: " << t << "    VCell: " << VCell << std::endl;
    //     }

    //     t++;

    // }
    

    // // // Multiply Grid Functions for Error Calculations
    // // CnCH_gf *= *domain_parameters.psi;
    // // phP_gf *= *domain_parameters.psi;

    // // CnE_gf *= *domain_parameters.pse;
    // // phE_gf *= *domain_parameters.pse;

    // // // Save simulation outputs
    // geometry.parallelMesh->Save("Results/pmesh");
    // domain_parameters.psi->Save("Results/psi");
    // domain_parameters.pse->Save("Results/pse");

    // // // // // CnP_gf.Save("Results/CnP");
    // CnCH_gf.Save("Results/CnCH");
    // CnE_gf.Save("Results/CnE");
    // // // phP_gf.Save("Results3/phP");
    // // // phE_gf.Save("Results3/phE");
    // Rxn_gf.Save("Results/Rxn");

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


