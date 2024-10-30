#include "CnP.hpp"
#include "mfem.hpp"


CnP::CnP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh) 
    
    {

    PsVc = mfem::HypreParVector(fespace); 


    }

void CnP::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation)
{
    Concentrations::SetInitialValues(Cn, initial_value, psx, perform_lithiation);

    Mmatp = std::make_shared<mfem::HypreParMatrix>();
    Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi);

    Concentrations::SetUpSolver(psx, Mmatp, *Mp_solver, Mp_prec);

    psx.GetTrueDofs(PsVc);
}
