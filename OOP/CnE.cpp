#include "CnE.hpp"
#include "mfem.hpp"


CnE::CnE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh) 
    
    {

    PeR = new ParGridFunction(fespace);
    RxE = new ParGridFunction(fespace);

    Mmate = std::make_shared<mfem::HypreParMatrix>();
    Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Me_prec.SetType(mfem::HypreSmoother::Jacobi);



    }

void CnE::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation)
{
    Concentrations::SetInitialValues(Cn, initial_value, psx, perform_lithiation);
    Concentrations::SetUpSolver(psx, Mmate, *Me_solver, Me_prec);

    ImposeNeumannBC(psx, *PeR);
}

void CnE::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{

    Concentrations::CreateReaction(Rx, *RxE, (-1.0 * Constants::t_minus));
    Concentrations::TotalReaction(*RxE, eCrnt);

    // Values used in Force Term Function
    mfem::ConstantCoefficient nbcCoef(infx); 
    mfem::GridFunctionCoefficient matCoef_R(PeR);
    mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    mfem::Array<int> nbc_w_bdr(pmesh->bdr_attributes.Max());
	nbc_w_bdr = 0; 
	nbc_w_bdr[0] = 1;

    Concentrations::ForceTerm(*RxE, ftE, nbc_w_bdr, m_nbcCoef, true); // true since applying BCs
    std::shared_ptr<GridFunctionCoefficient> cDe = Diffusivity(psx, Cn, false); // false using other equation





}
