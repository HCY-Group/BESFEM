#include "PotE.hpp"
#include "mfem.hpp"


PotE::PotE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh)
    
    {

    cgPE_solver = new mfem::CGSolver(MPI_COMM_WORLD);


    }

void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    Potentials::SetInitialPotentials(ph, initial_value);
    Potentials::SetUpSolver(*cgPE_solver, 1e-7, 80);

    Vcell = BvP - BvE;

    // std::cout << "Vcell: " << Vcell << std::endl;
    
}

