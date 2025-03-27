/**
 * @file CnP.cpp
 * @brief Implementation of the particle concentration class for battery simulations.
 */

#include "CnP.hpp"
#include "../code/Constants.hpp"
#include "mfem.hpp"
#include <optional>


using namespace std;

CnP::CnP(Initialize_Geometry &geo, Domain_Parameters &para)
    : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace),
    Tmatp(nullptr)
    
    {
    PsVc = mfem::HypreParVector(fespace.get());
    RxP = std::make_unique<mfem::ParGridFunction>(fespace.get());

    Mmatp = std::make_shared<mfem::HypreParMatrix>();
    Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi);

    Kmatp = std::make_shared<mfem::HypreParMatrix>();
    Fcb = mfem::HypreParVector(fespace.get());

    CpV0 = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    RHCp = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    CpVn = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));

    pKx2 = std::make_shared<mfem::ParBilinearForm>(fespace.get());

    }

void CnP::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    Concentrations::SetInitialConcentration(Cn, initial_value);
    Concentrations::LithiationCalculation(Cn, psx);
    Concentrations::SetUpSolver(psx, Mmatp, *Mp_solver, Mp_prec); // sets up Mass Matrix & Conditions

    psx.GetTrueDofs(PsVc); // Extract true degrees of freedom in the potential field

}

void CnP::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{
    // Compute the reaction field scaled by a constant factor
    Concentrations::CreateReaction(Rx, *RxP, (1.0/Constants::rho));

    // Assemble the force term without applying boundary conditions
    Solver::ForceTerm(*RxP, ftPC);
    
    // Compute the diffusivity coefficient and assemble the stiffness matrix
    std::shared_ptr<mfem::GridFunctionCoefficient> cDp = Concentrations::Diffusivity(psx, Cn, true); // true since using first equation
    pKx2 = std::make_shared<mfem::ParBilinearForm>(fespace.get());

    // Concentrations::KMatrix(boundary_dofs, Cn, ftPC, Kmatp, Fcb, cDp);
    Solver::StiffnessMatrix(cDp, boundary_dofs, Cn, ftPC, Kmatp, Fcb);
    pKx2->Update(fespace.get());

    Tmatp.reset(Add(1.0, *Mmatp, -(Constants::dt), *Kmatp));

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

    Cn.Distribute(CpVn.get());

    // Degree of Lithiation
    Concentrations::LithiationCalculation(Cn, psx);

}
