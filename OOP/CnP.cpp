#include "CnP.hpp"
#include "mfem.hpp"


CnP::CnP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh) 
    
    {

    PsVc = mfem::HypreParVector(fespace); 
    RxP = mfem::ParGridFunction(fespace);

    Mmatp = std::make_shared<mfem::HypreParMatrix>();
    Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi);



    }

void CnP::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation)
{
    Concentrations::SetInitialValues(Cn, initial_value, psx, perform_lithiation);
    Concentrations::SetUpSolver(psx, Mmatp, *Mp_solver, Mp_prec);

    psx.GetTrueDofs(PsVc);
}


void CnP::TimeStep(mfem::ParGridFunction &Rx)
{

    Concentrations::CreateReaction(Rx, RxP, (1/Constants::rho));






}
