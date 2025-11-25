#include "../include/CnC.hpp"
#include "../include/Constants.hpp"
#include "mfem.hpp"
#include <optional>


using namespace std;

CnC::CnC(Initialize_Geometry &geo, Domain_Parameters &para)
    : ConcentrationBase(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), fem(geo.parfespace), utils(geo,para), gtPsC(para.gtPsC), gtPsi(para.gtPsi),
    RxC(fespace.get()), Dp(fespace.get()), Mp_solver(MPI_COMM_WORLD), cAp(&RxC), cDp(&Dp), Fct(fespace.get()),
    PsVc(fespace.get()), CpV0(fespace.get()), RHCp(fespace.get()), CpVn(fespace.get())
    
    {

    if (gtPsC < 1.0e-200){
        gtPsC = gtPsi;
    }

    }

void CnC::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    utils.SetInitialValue(Cn, initial_value);
    utils.CalculateLithiation(Cn, psx, gtPsC);
    Xfr = utils.GetLithiation();


    mfem::GridFunctionCoefficient coef(&psx);
    fem.InitializeMassMatrix(coef, Mt);
    fem.FormSystemMatrix(Mt, boundary_dofs, Mmatp);

    Mp_solver.iterative_mode = false; // Enable iterative mode for the solver
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi); //
    fem.SolverConditions(Mmatp, Mp_solver, Mp_prec); // Set up the solver conditions for the mass matrix

    // fem.InitializeForceTerm(cAp, Bc2); // HALF
    // Fct = *Bc2; // Move the updated force term to Fct HALF

    fem.InitializeStiffnessMatrix(cDp, Kc2); // Initialize stiffness form for particle potential HALF

    psx.GetTrueDofs(PsVc); // Extract true degrees of freedom in the potential field

}

void CnC::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{
    // Compute the reaction field scaled by a constant factor
    utils.InitializeReaction(Rx, RxC, (1.0/Constants::rho_C));
    cAp.SetGridFunction(&RxC); // Set the reaction term coefficient for the force term

    // std::cout << "RxC Sum before: " << RxC.Sum() << std::endl;

    fem.InitializeForceTerm(cAp, Bc2);
    fem.Update(Bc2); // Update the force term with the current reaction term
    Fct = *Bc2; // Move the updated force term to Fct

    for (int vi = 0; vi < nV; vi++){
        Dp(vi) = psx(vi) * (0.0277 - 0.084 * Cn(vi) + 0.1003 * Cn(vi) * Cn(vi)) * 1.0e-8;
        if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}
    }
    cDp.SetGridFunction(&Dp); // Set the diffusivity coefficient for the stiffness matrix
    
    fem.Update(Kc2); // Update the stiffness matrix with the current diffusivity coefficient
    fem.FormLinearSystem(Kc2, boundary_dofs, Cn, Fct, Kmatp, X1v, Fcb); // Form the linear system for particle potential
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
    utils.CalculateLithiation(Cn, psx, gtPsC);
    Xfr = utils.GetLithiation();


}
