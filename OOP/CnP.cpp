/**
 * @file CnP.cpp
 * @brief Implementation of the particle concentration class for battery simulations.
 */

#include "CnP.hpp"
#include "mfem.hpp"

CnP::CnP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh), fespace(fe)
    
    {
    PsVc = mfem::HypreParVector(fespace);
    RxP = new mfem::ParGridFunction(fespace);

    Mmatp = std::make_shared<mfem::HypreParMatrix>();
    Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi);

    Kmatp = std::make_shared<mfem::HypreParMatrix>();
    Fcb = HypreParVector(fespace);

    CpV0 = new mfem::HypreParVector(fespace);
    RHCp = new mfem::HypreParVector(fespace);
    CpVn = new mfem::HypreParVector(fespace);

    // std::cout << "fespace address in CnP: " << fespace << std::endl;

    pKx2 = std::make_shared<mfem::ParBilinearForm>(fespace);
    // std::cout << "fespace in CnP Constr: " << fespace << std::endl;

    

    }

void CnP::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    Concentrations::SetInitialConcentration(Cn, initial_value);
    Concentrations::LithiationCalculation(Cn, psx);
    Concentrations::SetUpSolver(psx, Mmatp, *Mp_solver, Mp_prec);

    psx.GetTrueDofs(PsVc); // Extract true degrees of freedom in the potential field

    // pKx2 = std::make_shared<mfem::ParBilinearForm>(fespace);

}


void CnP::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{
    // Compute the reaction field scaled by a constant factor
    Concentrations::CreateReaction(Rx, *RxP, (1.0/Constants::rho));

    // Create dummy values to use Force Term Function 
    static Array<int> dummy_boundary;
    static mfem::ConstantCoefficient coef1(0.0);
    static mfem::ConstantCoefficient coef2(0.0);
    static mfem::ProductCoefficient dummy_coef(coef1, coef2);

    // Assemble the force term without applying boundary conditions
    Concentrations::ForceTerm(*RxP, ftPC, dummy_boundary, dummy_coef, false); // false since not applying BCs

    pKx2->Update(fespace);

    // Compute the diffusivity coefficient and assemble the stiffness matrix
    std::shared_ptr<GridFunctionCoefficient> cDp = Concentrations::Diffusivity(psx, Cn, true); // true since using first equation
    pKx2 = std::make_shared<mfem::ParBilinearForm>(fespace);


    Concentrations::KMatrix(pKx2, boundary_dofs, Cn, ftPC, Kmatp, X1v, Fcb, cDp);

    // cDp->ResetCoefficient();
    pKx2->Update(fespace);


    // Form the time-stepping system matrix
    Tmatp = Add(1.0, *Mmatp, -(Constants::dt), *Kmatp);

    // Solve for the next time-step concentrations
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
