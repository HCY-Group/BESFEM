

#include "../include/CnC.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <optional>


using namespace std;

CnC::CnC(Initialize_Geometry &geo, Domain_Parameters &para)
    : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), gtPsC(para.gtPsC), gtPsi(para.gtPsi),
    RxC(fespace.get()), Dp(fespace.get()), Mp_solver(MPI_COMM_WORLD), Fct(fespace.get()), cAp(&RxC), cDp(&Dp),
    PsVc(fespace.get()), CpV0(fespace.get()), RHCp(fespace.get()), CpVn(fespace.get())
    
    {

    // std::cout << "gtPsC before: " << gtPsC << std::endl;

    if (gtPsC < 1.0e-200){
        gtPsC = gtPsi;
    }

    }

void CnC::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    Concentrations::SetInitialConcentration(Cn, initial_value);
    Concentrations::LithiationCalculation(Cn, psx, gtPsC);

    mfem::GridFunctionCoefficient coef(&psx);
    SolverSteps::InitializeMassMatrix(coef, Mt);
    SolverSteps::FormSystemMatrix(Mt, boundary_dofs, Mmatp);

    Mp_solver.iterative_mode = false; // Enable iterative mode for the solver
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi); //
    SolverSteps::SolverConditions(Mmatp, Mp_solver, Mp_prec); // Set up the solver conditions for the mass matrix

    // SolverSteps::InitializeForceTerm(cAp, Bc2); // HALF
    // Fct = *Bc2; // Move the updated force term to Fct HALF

    SolverSteps::InitializeStiffnessMatrix(cDp, Kc2); // Initialize stiffness form for particle potential HALF

    psx.GetTrueDofs(PsVc); // Extract true degrees of freedom in the potential field

}

void CnC::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{
    // Compute the reaction field scaled by a constant factor
    Concentrations::CreateReaction(Rx, RxC, (1.0/Constants::rho_C));
    cAp.SetGridFunction(&RxC); // Set the reaction term coefficient for the force term

    SolverSteps::InitializeForceTerm(cAp, Bc2);
    SolverSteps::Update(Bc2); // Update the force term with the current reaction term
    Fct = *Bc2; // Move the updated force term to Fct

    for (int vi = 0; vi < nV; vi++){
        Dp(vi) = psx(vi) * (0.0277 - 0.084 * Cn(vi) + 0.1003 * Cn(vi) * Cn(vi)) * 1.0e-8;
        if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}
    }
    cDp.SetGridFunction(&Dp); // Set the diffusivity coefficient for the stiffness matrix
    
    SolverSteps::Update(Kc2); // Update the stiffness matrix with the current diffusivity coefficient
    SolverSteps::FormLinearSystem(Kc2, boundary_dofs, Cn, Fct, Kmatp, X1v, Fcb); // Form the linear system for particle potential
    Fcb *= Constants::dt; // Scale the right-hand side vector by the time step

    Tmatp.reset(Add(1.0, Mmatp, -1.0* Constants::dt, Kmatp));

    // Solve for the next time-step concentrations
    Cn.GetTrueDofs(CpV0);
    Tmatp->Mult(CpV0, RHCp);
    RHCp += Fcb;

    Mp_solver.Mult(RHCp, CpVn);

    // Update only the solid region MAKE INTO FUNCTION
    for (int p = 0; p < CpV0.Size(); p++){
        if (PsVc(p) < 1.0e-5){
            (CpVn)(p) = Constants::init_CnC;} // Cp0 initial value
    }

    Cn.Distribute(CpVn);

    // Degree of Lithiation
    Concentrations::LithiationCalculation(Cn, psx, gtPsC);

}
