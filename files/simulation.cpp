#include "mfem.hpp"
#include "mpi.h"

#include "Constants.hpp"
#include "Mesh_Handler.hpp"
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

    using namespace std::chrono;
    auto program_start = high_resolution_clock::now();


    // Initialize MPI and HYPRE
    Mpi::Init(argc, argv);
    Hypre::Init();

    // Create the MeshHandler object
    MeshHandler mesh_handler;
    mesh_handler.LoadMesh();

    // Create pmesh & fespace to use for all 
    ParMesh pmesh = mesh_handler.GetMesh();
    H1_FECollection fec(Constants::order, pmesh.Dimension());
    ParFiniteElementSpace fespace(&pmesh, &fec);
    mesh_handler.SetupBoundaryConditions(&pmesh, &fespace);

    // Retrieve psi & pse from mesh_handler
    mfem::ParGridFunction &psi = *mesh_handler.GetPsi();
    mfem::ParGridFunction &pse = *mesh_handler.GetPse();

    // Initialize CnP & CnE
    CnP particle_concentration(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction CnP_gf(&fespace);
    particle_concentration.Initialize(CnP_gf, 0.3, psi); 

    CnE electrolyte_concentration(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction CnE_gf(&fespace);
    electrolyte_concentration.Initialize(CnE_gf, 0.001, pse); 

    // Initialize phP & phE
    PotP particle_potential(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction phP_gf(&fespace);
    particle_potential.Initialize(phP_gf, 2.9395);

    PotE electrolyte_potential(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction phE_gf(&fespace);
    electrolyte_potential.Initialize(phE_gf, -1.0);
    
    // Initialize Reaction
    Reaction reaction(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction Rxn_gf(&fespace);
    reaction.Initialize(Rxn_gf, 0.0);

    // Create Current Class
    Current current(&pmesh, &fespace, mesh_handler);

    double global_current = 0.0;
    double VCell = BvP - BvE;
    // std::cout << "VCell Initial: " << VCell << std::endl;

    // std::cout << "VCut: " << Constants::VCut << std::endl;
 
    // Time Step
    for (int t = 0; t < 20 + 1; ++t) {
    // while ( VCell > Constants::VCut) {

        particle_concentration.TimeStep(Rxn_gf, CnP_gf, psi);
        electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, pse);
        
        particle_potential.TimeStep(CnP_gf, psi, phP_gf);
        electrolyte_potential.TimeStep(CnE_gf, pse, phE_gf);

        // rate constants and exchange current density at interface
        reaction.ExchangeCurrentDensity(CnP_gf); 

        // while loop
        double globalerror_P = 1.0;
		double globalerror_E = 1.0;

        // int inlp = 0;
        while (globalerror_P > 1.0e-9 || globalerror_E > 1.0e-9) {
            reaction.ButlerVolmer(Rxn_gf, CnP_gf, CnE_gf, phP_gf, phE_gf);
            particle_potential.CalculateGlobalError(Rxn_gf, phP_gf, psi, globalerror_P);
            electrolyte_potential.CalculateGlobalError(Rxn_gf, phE_gf, pse, globalerror_E);
        }

        // total reaction current
        reaction.TotalReactionCurrent(Rxn_gf, global_current);
        current.Constant(phP_gf, global_current);

        VCell = BvP - BvE;
        // std::cout << "VCell: " << VCell << std::endl;


    }

    // pmesh.Save("pmesh");
    // CnP_gf.Save("CnP");
    // CnE_gf.Save("CnE");
    // phP_gf.Save("phP");
    // phE_gf.Save("phE");

    Mpi::Finalize();

    // End timing the entire program
    auto program_end = high_resolution_clock::now();
    std::cout << "Total Program Time: " 
              << duration_cast<seconds>(program_end - program_start).count() 
              << " seconds" << std::endl;



    return 0;
}


