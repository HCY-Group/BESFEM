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

// #include "Reaction.hpp"
// #include "Potentials.hpp"

#include <iostream>

int main(int argc, char *argv[]) {
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
    particle_concentration.Initialize(CnP_gf, 0.3, psi, true); // true to run lithiation calculation

    CnE electrolyte_concentration(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction CnE_gf(&fespace);
    electrolyte_concentration.Initialize(CnE_gf, 0.001, pse, false); // false since not running lithiation calculation

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
 
    // Time Step
    for (int t = 0; t < 10 + 1; ++t) {

        particle_concentration.TimeStep(Rxn_gf, CnP_gf, psi);
        electrolyte_concentration.TimeStep(Rxn_gf, CnE_gf, pse);
        
        particle_potential.TimeStep(CnP_gf, psi, phP_gf);
        electrolyte_potential.TimeStep(CnE_gf, pse, phE_gf);

        // rate constants and exchange current density at interface
        reaction.ExchangeCurrentDensity(CnP_gf); // move to CnP? kfw and kbw used in reaction BV

        // while loop
        reaction.ButlerVolmer(Rxn_gf, CnP_gf, CnE_gf, phP_gf, phE_gf);
        // particle_potential.CalculateGlobalError(Rxn_gf, phP_gf, psi);
        // std::cout << "Rxn: " << Rxn_gf << std::endl;
        // Rxn_gf.Print(std::cout);

    }

    // particle_concentration.Save(CnP_gf, "CnP");
    // electrolyte_concentration.Save(CnE_gf, "CnE");

    Mpi::Finalize();
    return 0;
}


