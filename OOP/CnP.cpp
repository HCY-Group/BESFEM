#include "CnP.hpp"
#include "mfem.hpp"


CnP::CnP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh)
    
    {

    PsVc = mfem::HypreParVector(fespace); 
    RxP = new ParGridFunction(fespace);

    Mmatp = std::make_shared<mfem::HypreParMatrix>();
    Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi);

    Kmatp = std::make_shared<mfem::HypreParMatrix>();
    Fcb = HypreParVector(fespace);

    CpV0 = new mfem::HypreParVector(fespace);
    RHCp = new mfem::HypreParVector(fespace);
    CpVn = new mfem::HypreParVector(fespace);


    }

void CnP::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation)
{
    Concentrations::SetInitialValues(Cn, initial_value, psx, perform_lithiation);
    Concentrations::SetUpSolver(psx, Mmatp, *Mp_solver, Mp_prec);

    psx.GetTrueDofs(PsVc);
}


void CnP::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{

    Concentrations::CreateReaction(Rx, *RxP, (1/Constants::rho));

    // Create dummy values to use Force Term Function 
    Array<int> dummy_boundary;
    mfem::ConstantCoefficient coef1(0.0);
    mfem::ConstantCoefficient coef2(0.0);
    mfem::ProductCoefficient dummy_coef(coef1, coef2);

    Concentrations::ForceTerm(*RxP, ftP, dummy_boundary, dummy_coef, false); // false since not applying BCs

    std::shared_ptr<GridFunctionCoefficient> cDp = Concentrations::Diffusivity(psx, Cn, true); // true since using first equation
    Concentrations::KMatrix(boundary_dofs, Cn, ftP, Kmatp, X1v, Fcb, cDp.get());

    Tmatp = Add(1.0, *Mmatp, -(Constants::dt), *Kmatp);

    int nDof = CpV0->Size();
    Cn.GetTrueDofs(*CpV0);

    Tmatp->Mult(*CpV0, *RHCp);
    *RHCp += Fcb;

    Mp_solver->Mult(*RHCp, *CpVn);

    // Update only the solid region MAKE INTO FUNCTION
    for (int p = 0; p < nDof; p++){
        if (PsVc(p) < 1.0e-5){
            (*CpVn)(p) = 0.3;} // Cp0 initial value
    }

    Cn.Distribute(CpVn);


    // Degree of Lithiation
    Concentrations::LithiationCalculation(Cn, psx);


}

void CnP::Save(mfem::ParGridFunction &gf, const std::string &base_name){

    Concentrations::Save(gf, base_name);

}
