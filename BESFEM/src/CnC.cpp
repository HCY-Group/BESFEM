

#include "../include/CnC.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <optional>


using namespace std;

CnC::CnC(Initialize_Geometry &geo, Domain_Parameters &para)
    : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace),
    RxP(fespace.get()), Dp(fespace.get()), Mp_solver(MPI_COMM_WORLD), Fct(fespace.get()), cAp(&RxP), cDp(&Dp)
    
    {
    PsVc = mfem::HypreParVector(fespace.get());

    Mmatp = std::make_shared<mfem::HypreParMatrix>();
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi);

    Kmatp = std::make_shared<mfem::HypreParMatrix>();

    CpV0 = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    RHCp = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    CpVn = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));

    }

void CnC::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    Concentrations::SetInitialConcentration(Cn, initial_value);
    Concentrations::LithiationCalculation(Cn, psx);

    mfem::GridFunctionCoefficient coef(&psx);
    SolverSteps::InitializeMassMatrix(coef, Mt);
    SolverSteps::FormSystemMatrix(Mt, boundary_dofs, *Mmatp);

    Mp_solver.iterative_mode = false; // Enable iterative mode for the solver
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi); //
    SolverSteps::SolverConditions(*Mmatp, Mp_solver, Mp_prec); // Set up the solver conditions for the mass matrix

    SolverSteps::InitializeForceTerm(cAp, Bc2);
    Fct = *Bc2; // Move the updated force term to Fct

    SolverSteps::InitializeStiffnessMatrix(cDp, Kc2); // Initialize stiffness form for particle potential

    psx.GetTrueDofs(PsVc); // Extract true degrees of freedom in the potential field

}

void CnC::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{
    // Compute the reaction field scaled by a constant factor
    Concentrations::CreateReaction(Rx, RxP, (1.0/Constants::rho_C));
    cAp.SetGridFunction(&RxP); // Set the reaction term coefficient for the force term

    // SolverSteps::InitializeForceTerm(cAp, Bc2);
    SolverSteps::Update(Bc2); // Update the force term with the current reaction term
    Fct = *Bc2; // Move the updated force term to Fct

    // Compute the diffusivity coefficient and assemble the stiffness matrix
    std::shared_ptr<mfem::GridFunctionCoefficient> cDp = Concentrations::Diffusivity(psx, Cn, true); // true since using first equation

    SolverSteps::Update(Kc2); // Update the stiffness matrix with the current diffusivity coefficient
    SolverSteps::FormLinearSystem(Kc2, boundary_dofs, Cn, Fct, *Kmatp, X1v, Fcb); // Form the linear system for particle potential
    Fcb *= Constants::dt; // Scale the right-hand side vector by the time step

    Tmatp.reset(Add(1.0, *Mmatp, -(Constants::dt), *Kmatp));

    // Solve for the next time-step concentrations
    int nDof = CpV0->Size();
    Cn.GetTrueDofs(*CpV0);
    Tmatp->Mult(*CpV0, *RHCp);
    *RHCp += Fcb;

    Mp_solver.Mult(*RHCp, *CpVn);

    // Update only the solid region MAKE INTO FUNCTION
    for (int p = 0; p < nDof; p++){
        if (PsVc(p) < 1.0e-5){
            (*CpVn)(p) = 0.3;} // Cp0 initial value
    }

    Cn.Distribute(CpVn.get());

    // Degree of Lithiation
    Concentrations::LithiationCalculation(Cn, psx);

}
