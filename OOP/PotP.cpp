#include "PotP.hpp"
#include "mfem.hpp"


PotP::PotP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh)
    
    {

    cgPP_solver = new mfem::CGSolver(MPI_COMM_WORLD);


    }

void PotP::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    Potentials::SetInitialPotentials(ph, initial_value);
    Potentials::SetUpSolver(*cgPP_solver, 1e-7, 82);

    
}

