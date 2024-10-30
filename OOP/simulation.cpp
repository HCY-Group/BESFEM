#include "mfem.hpp"
#include "mpi.h"

#include "Constants.hpp"
#include "Mesh_Handler.hpp"
#include "Concentrations_Base.hpp"
#include "CnP.hpp"
#include "CnE.hpp"

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

    // Retrieve psi from mesh_handler
    mfem::ParGridFunction &psi = *mesh_handler.GetPsi();

    // Define Mmatp, Mp_solver, and Mp_prec for particle concentration
    std::shared_ptr<mfem::HypreParMatrix> Mmatp = std::make_shared<mfem::HypreParMatrix>();
    std::shared_ptr<mfem::CGSolver> Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    mfem::HypreSmoother Mp_prec; 

    // Retrieve pse from mesh_handler
    mfem::ParGridFunction &pse = *mesh_handler.GetPse();

    // Define Mmate, Me_solver, and Me_prec for electrolyte concentration
    std::shared_ptr<mfem::HypreParMatrix> Mmate = std::make_shared<mfem::HypreParMatrix>();
    std::shared_ptr<mfem::CGSolver> Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    mfem::HypreSmoother Me_prec;

    // Initialize CnP & CnE
    CnP particle_concentration(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction CnP_gf(&fespace);
    particle_concentration.Initialize(CnP_gf, 0.3, psi, Mmatp, *Mp_solver, Mp_prec, true); // true to run lithiation 

    CnE electrolyte_concentration(&pmesh, &fespace, mesh_handler);
    mfem::ParGridFunction CnE_gf(&fespace);
    electrolyte_concentration.Initialize(CnE_gf, 0.001, pse, Mmate, *Me_solver, Me_prec, false);

    Mpi::Finalize();
    return 0;
}


