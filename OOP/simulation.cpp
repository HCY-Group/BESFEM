#include "mfem.hpp"
#include "mpi.h"

#include "Constants.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "Concentrations_Base.hpp"
#include "Potentials_Base.hpp"
#include "CnP.hpp"
#include "CnE.hpp"
#include "PotP.hpp"
#include "PotE.hpp"
#include "Reaction.hpp"
#include "Current.hpp"

#include <chrono>
#include <iostream>

int main(int argc, char *argv[]) {

    // Start measuring the program execution time
    using namespace std::chrono;
    auto program_start = high_resolution_clock::now();

    // Initialize MPI for parallel processing and HYPRE for solver setup
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    // // Create the MeshHandler object and load the mesh
    // MeshHandler mesh_handler;
    // mesh_handler.LoadMesh();

    // Initialize Mesh & Geometry
    Initialize_Geometry geometry;
    geometry.InitializeMesh(Constants::mesh_file, MPI_COMM_WORLD, Constants::order);
    geometry.SetupBoundaryConditions();

    // Initialize and Calculate Domain Parameters (psi, pse, AvB, AvP)
    Domain_Parameters domain_parameters(geometry);
    domain_parameters.SetupDomainParameters();

    // // Create a parallel mesh and finite element space to use across the simulation 
    // mfem::ParMesh &pmesh = *mesh_handler.pmesh.get();
    // mfem::H1_FECollection fec(Constants::order, pmesh.Dimension());
    // ParFiniteElementSpace fespace(&pmesh, &fec);

    // // Setup boundary conditions for the mesh using the MeshHandler
    // mesh_handler.SetupBoundaryConditions(&pmesh, &fespace);

    // // Retrieve psi (particle phase potential) and pse (electrolyte phase potential) from the MeshHandler
    // mfem::ParGridFunction &psi = *mesh_handler.psi;
    // mfem::ParGridFunction &pse = *mesh_handler.pse;

    // Initialize Concentration Classes

    // Particle concentration initialization with a grid function and initial value
    CnP particle_concentration(geometry, domain_parameters);
    mfem::ParGridFunction CnP_gf(geometry.parfespace.get());
    particle_concentration.Initialize(CnP_gf, 0.3, *domain_parameters.psi); 

    // Electrolyte concentration initialization with a grid function and initial value
    CnE electrolyte_concentration(geometry, domain_parameters);
    mfem::ParGridFunction CnE_gf(geometry.parfespace.get());
    electrolyte_concentration.Initialize(CnE_gf, 0.001, *domain_parameters.pse); 

    // Initialize Potential Classes

    // Particle potential initialization with a grid function and initial value
    PotP particle_potential(geometry, domain_parameters);
    mfem::ParGridFunction phP_gf(geometry.parfespace.get());
    particle_potential.Initialize(phP_gf, 2.9395);

    // Electrolyte potential initialization with a grid function and initial value
    PotE electrolyte_potential(geometry, domain_parameters);
    mfem::ParGridFunction phE_gf(geometry.parfespace.get());
    electrolyte_potential.Initialize(phE_gf, -1.0);
    
    // Initialize Reaction Class

    // Initialize the reaction with an empty reaction grid function and initial value
    Reaction reaction(geometry, domain_parameters);
    mfem::ParGridFunction Rxn_gf(geometry.parfespace.get());
    reaction.Initialize(Rxn_gf, 0.0);

    // Initialize Current Class

    // Create the Current class to control current based on particle potential
    Current current(geometry, domain_parameters);

    // Set initial global current and cell voltage values  
    double global_current = 0.0;
    double VCell = BvP - BvE;

    // Main Simulation Loop
 
    // Perform simulation over time steps
    for (int t = 0; t < 6000 + 1; ++t) {
    // while ( VCell > Constants::VCut) {

        // Step 1: Update concentrations for both particle and electrolyte phases
        particle_concentration.TimeStep(Rxn_gf, CnP_gf, *domain_parameters.psi);
        electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, *domain_parameters.pse);
        
        // Step 2: Update potentials for both particle and electrolyte phases
        particle_potential.TimeStep(CnP_gf, *domain_parameters.psi, phP_gf);
        electrolyte_potential.TimeStep(CnE_gf, *domain_parameters.pse, phE_gf);

        // Step 3: Compute rate constants and exchange current density at the interface     
        reaction.ExchangeCurrentDensity(CnP_gf); 

        // Step 4: Iteratively solve for global reaction rates using Butler-Volmer kinetics
        double globalerror_P = 1.0; // Error for particle potential
		double globalerror_E = 1.0; // Error for electrolyte potential

        // int inlp = 0;
        while (globalerror_P > 1.0e-9 || globalerror_E > 1.0e-9) {

            // Update reaction rates using the Butler-Volmer equation
            reaction.ButlerVolmer(Rxn_gf, CnP_gf, CnE_gf, phP_gf, phE_gf);

            // Calculate global errors for particle and electrolyte potentials
            particle_potential.CalculateGlobalError(Rxn_gf, phP_gf, *domain_parameters.psi, globalerror_P);
            electrolyte_potential.CalculateGlobalError(Rxn_gf, phE_gf, *domain_parameters.pse, globalerror_E);
        }

        // Step 5: Compute the total reaction current
        reaction.TotalReactionCurrent(Rxn_gf, global_current);

        // Step 6: Adjust the particle potential based on the target constant current
        current.Constant(phP_gf, global_current);

        // Step 7: Update the cell voltage based on the particle and electrolyte potentials
        VCell = BvP - BvE;
        
        std::cout << "VCell: " << VCell << std::endl;

        // if (t % 100 == 0) {
        // std::cout << "Iteration " << t << ": VCell = " << VCell << std::endl;
        // }

    }

    // Save simulation outputs
    // pmesh.Save("Results/pmesh");
    CnP_gf.Save("Results/CnP");
    CnE_gf.Save("Results/CnE");
    phP_gf.Save("Results/phP");
    phE_gf.Save("Results/phE");
    // psi.Save("Results/psi");
    // pse.Save("Results/pse");

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    // End timing and output the total program execution time
    auto program_end = high_resolution_clock::now();
    std::cout << "Total Program Time: " 
              << duration_cast<seconds>(program_end - program_start).count() 
              << " seconds" << std::endl;



    return 0;
}


