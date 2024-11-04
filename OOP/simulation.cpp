#include "mfem.hpp"
#include "mpi.h"

#include "Constants.hpp"
#include "Mesh_Handler.hpp"
#include "Concentrations_Base.hpp"
#include "CnP.hpp"
#include "CnE.hpp"
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

    // Initialize Reaction
    Reaction reaction(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction Rxn_gf(&fespace);
    reaction.Initialize(Rxn_gf, 0.0);
 
    // Time Step
    for (int t = 0; t < 10 + 1; ++t) {
        particle_concentration.TimeStep(Rxn_gf);
        electrolyte_concentration.TimeStep(Rxn_gf);


    }

    Mpi::Finalize();
    return 0;
}


