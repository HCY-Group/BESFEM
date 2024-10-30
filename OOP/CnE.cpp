#include "CnE.hpp"
#include "mfem.hpp"


CnE::CnE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh) 
    
    {

    PeR = new ParGridFunction(fespace);



    }

void CnE::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation)
{
    Concentrations::SetInitialValues(Cn, initial_value, psx, perform_lithiation);

    Mmate = std::make_shared<mfem::HypreParMatrix>();
    Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Me_prec.SetType(mfem::HypreSmoother::Jacobi);

    Concentrations::SetUpSolver(psx, Mmate, *Me_solver, Me_prec);

    ImposeNeumannBC(psx, *PeR);
}
