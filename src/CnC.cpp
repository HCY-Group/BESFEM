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

// void CnC::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
// {
//     utils.SetInitialValue(Cn, initial_value);
//     utils.CalculateLithiation(Cn, psx, gtPsC);
//     Xfr = utils.GetLithiation();


//     mfem::GridFunctionCoefficient coef(&psx);
//     fem.InitializeMassMatrix(coef, Mt);
//     fem.FormSystemMatrix(Mt, boundary_dofs, Mmatp);

//     Mp_solver.iterative_mode = false; // Enable iterative mode for the solver
//     Mp_prec.SetType(mfem::HypreSmoother::Jacobi); //
//     fem.SolverConditions(Mmatp, Mp_solver, Mp_prec); // Set up the solver conditions for the mass matrix

//     // fem.InitializeForceTerm(cAp, Bc2); // HALF
//     // Fct = *Bc2; // Move the updated force term to Fct HALF

//     fem.InitializeStiffnessMatrix(cDp, Kc2); // Initialize stiffness form for particle potential HALF

//     psx.GetTrueDofs(PsVc); // Extract true degrees of freedom

// }

void CnC::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx)
{
    boundary_dofs.DeleteAll();

    utils.SetInitialValue(Cn, initial_value);
    utils.CalculateLithiation(Cn, psx, gtPsx);
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

    psx.GetTrueDofs(PsVc); // Extract true degrees of freedom

    // // testing purposes
    // utils.CalculateLithiation(Cn, psx, gtPsx);
    // Xfr = utils.GetLithiation();

    // std::cout << "Initial Lithiation Fraction: " << Xfr << std::endl;


}

// void CnC::Particle_Particle(mfem::ParGridFunction &sum_part, mfem::ParGridFunction &weight, mfem::ParGridFunction &grad_psi, mfem::ParGridFunction &mu_1, mfem::ParGridFunction &mu_2)
// {

//     for (int vi = 0; vi < nV; vi++){

//         double grad_psi_val = grad_psi(vi);
//         double weight_val = weight(vi);
//         double mu1_val = mu_1(vi);
//         double mu2_val = mu_2(vi);

//         sum_part(vi) = weight_val * grad_psi_val * Constants::rho_C * (1.0/Constants::RT) * Constants::Perm * (mu2_val - mu1_val);
//     }

//     // sum_part.SaveAsOne("sum_part");

// }


// void CnC::UpdateConcentration(mfem::ParGridFunction &Rx,mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtPsx, 
//                             mfem::ParGridFunction &weight_elec, 
//                             mfem::ParGridFunction &sum_AB, mfem::ParGridFunction &weight_AB, mfem::ParGridFunction &grad_AB,
//                             mfem::ParGridFunction &sum_AC, mfem::ParGridFunction &weight_AC, mfem::ParGridFunction &grad_AC,
//                             mfem::ParGridFunction &mu_A, mfem::ParGridFunction &mu_B, mfem::ParGridFunction &mu_C, mfem::ParGridFunction &mu_D, mfem::ParGridFunction &psiA, mfem::ParGridFunction &psiB, mfem::ParGridFunction &psiC)
// {
    
//     // Compute the reaction field scaled by a constant factor
//     // utils.InitializeReaction(Rx, RxC, (1.0/Constants::rho_C));

//     utils.InitializeReaction(Rx, RxC, (1.0));

//     RxC *= weight_elec; // Scale reaction by electrode weighting function

//     Particle_Particle(sum_AB, weight_AB, grad_AB, mu_A, mu_B);
//     Particle_Particle(sum_AC, weight_AC, grad_AC, mu_C, mu_D);

    


//     // sum_AB *= psiA;
//     // sum_AB *= psiB;

//     // sum_AC *= psiA;
//     // sum_AC *= psiC;

    
//     RxC += sum_AB;
//     RxC += sum_AC;

//     // RxC.SaveAsOne("RxC_total");



//     cAp.SetGridFunction(&RxC); // Set the reaction term coefficient for the force term

//     // RxC.SaveAsOne("RxC");

//     fem.InitializeForceTerm(cAp, Bc2);
//     fem.Update(Bc2); // Update the force term with the current reaction term
//     Fct = *Bc2; // Move the updated force term to Fct

//     for (int vi = 0; vi < nV; vi++){
//         Dp(vi) = psx(vi) * (0.0277 - 0.084 * Cn(vi) + 0.1003 * Cn(vi) * Cn(vi)) * 1.0e-8;
//         if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}
//         // std::cout << "Diffusivity at node " << vi << ": " << Dp(vi) << std::endl;
//     }
//     cDp.SetGridFunction(&Dp); // Set the diffusivity coefficient for the stiffness matrix
    
//     fem.Update(Kc2); // Update the stiffness matrix with the current diffusivity coefficient
//     fem.FormLinearSystem(Kc2, boundary_dofs, Cn, Fct, Kmatp, X1v, Fcb); // Form the linear system for particle potential
//     Fcb *= Constants::dt; // Scale the right-hand side vector by the time step

//     Tmatp.reset(Add(1.0, Mmatp, -1.0* Constants::dt, Kmatp));

//     // Solve for the next time-step concentrations
//     Cn.GetTrueDofs(CpV0);
//     Tmatp->Mult(CpV0, RHCp);
//     RHCp += Fcb;

//     Mp_solver.Mult(RHCp, CpVn);

//     psx.GetTrueDofs(PsVc);

//     mfem::Vector ps_true;
//     psx.GetTrueDofs(ps_true);

//     // Update only the solid region MAKE INTO FUNCTION
//     for (int p = 0; p < CpV0.Size(); p++){
//         if (PsVc(p) < 1.0e-5){ // 1e-1 works, but gaps still get smaller
//             // (CpVn)(p) = Constants::init_CnC;} // Cp0 initial value
//             (CpVn)(p) = CpV0(p);} // Cp0 initial value

//         if (CpVn(p) < 0.0){
//             (CpVn)(p) = 1.0e-6;} // prevent negative concentrations
//     }

//     Cn.Distribute(CpVn);

//     // Degree of Lithiation
//     utils.CalculateLithiation(Cn, psx, gtPsx);
//     Xfr = utils.GetLithiation();

//     // std::cout << "Updated Lithiation Fraction: " << Xfr << std::endl;

//     Rx = RxC;


// }

void CnC::ComputePairFlux(mfem::ParGridFunction &sum_part, mfem::ParGridFunction &weight, mfem::ParGridFunction &grad_psi, mfem::ParGridFunction &mu_1, mfem::ParGridFunction &mu_2)
{
    for (int vi = 0; vi < nV; vi++){

        double grad_psi_val = grad_psi(vi);
        double weight_val = weight(vi);
        double mu1_val = mu_1(vi);
        double mu2_val = mu_2(vi);

        sum_part(vi) = weight_val * grad_psi_val * Constants::rho_C * (1.0/Constants::RT) * Constants::Perm * (mu2_val - mu1_val);
    }

}

void CnC::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                            double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms)

{
    
    utils.InitializeReaction(Rx, RxC, (1.0));
    RxC *= weight_elec; // Scale reaction by electrode weighting function

    // Add all particle-particle pair exchange terms
    for (const auto &pair : pair_terms)
    {
        MFEM_VERIFY(pair.sum_part, "pair.sum_part is null");
        MFEM_VERIFY(pair.weight,   "pair.weight is null");
        MFEM_VERIFY(pair.grad_psi, "pair.grad_psi is null");
        MFEM_VERIFY(pair.mu_self,  "pair.mu_self is null");
        MFEM_VERIFY(pair.mu_nbr,   "pair.mu_nbr is null");

        ComputePairFlux(*pair.sum_part, *pair.weight, *pair.grad_psi, *pair.mu_self, *pair.mu_nbr);

        RxC += *pair.sum_part;
    }



    // Particle_Particle(sum_AB, weight_AB, grad_AB, mu_A, mu_B);
    // Particle_Particle(sum_AC, weight_AC, grad_AC, mu_C, mu_D);

    


    // sum_AB *= psiA;
    // sum_AB *= psiB;

    // sum_AC *= psiA;
    // sum_AC *= psiC;

    
    // RxC += sum_AB;
    // RxC += sum_AC;

    // RxC.SaveAsOne("RxC_total");



    cAp.SetGridFunction(&RxC); // Set the reaction term coefficient for the force term

    // RxC.SaveAsOne("RxC");

    fem.InitializeForceTerm(cAp, Bc2);
    fem.Update(Bc2); // Update the force term with the current reaction term
    Fct = *Bc2; // Move the updated force term to Fct

    for (int vi = 0; vi < nV; vi++){
        Dp(vi) = psx(vi) * (0.0277 - 0.084 * Cn(vi) + 0.1003 * Cn(vi) * Cn(vi)) * 1.0e-8;
        if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}
        // std::cout << "Diffusivity at node " << vi << ": " << Dp(vi) << std::endl;
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

    psx.GetTrueDofs(PsVc);

    mfem::Vector ps_true;
    psx.GetTrueDofs(ps_true);

    // Update only the solid region MAKE INTO FUNCTION
    for (int p = 0; p < CpV0.Size(); p++){
        if (PsVc(p) < 1.0e-5){ // 1e-1 works, but gaps still get smaller
            // (CpVn)(p) = Constants::init_CnC;} // Cp0 initial value
            (CpVn)(p) = CpV0(p);} // Cp0 initial value

        if (CpVn(p) < 0.0){
            (CpVn)(p) = 1.0e-6;} // prevent negative concentrations
    }

    Cn.Distribute(CpVn);

    // Degree of Lithiation
    utils.CalculateLithiation(Cn, psx, gtPsx);
    Xfr = utils.GetLithiation();

    // std::cout << "Updated Lithiation Fraction: " << Xfr << std::endl;

    Rx = RxC;


}

// void CnC::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
// {
//     // Compute the reaction field scaled by a constant factor
//     utils.InitializeReaction(Rx, RxC, (1.0/Constants::rho_C));
//     cAp.SetGridFunction(&RxC); // Set the reaction term coefficient for the force term

//     // std::cout << "RxC Sum before: " << RxC.Sum() << std::endl;

//     fem.InitializeForceTerm(cAp, Bc2);
//     fem.Update(Bc2); // Update the force term with the current reaction term
//     Fct = *Bc2; // Move the updated force term to Fct

//     for (int vi = 0; vi < nV; vi++){
//         Dp(vi) = psx(vi) * (0.0277 - 0.084 * Cn(vi) + 0.1003 * Cn(vi) * Cn(vi)) * 1.0e-8;
//         // if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}
//     }
//     cDp.SetGridFunction(&Dp); // Set the diffusivity coefficient for the stiffness matrix
    
//     fem.Update(Kc2); // Update the stiffness matrix with the current diffusivity coefficient
//     fem.FormLinearSystem(Kc2, boundary_dofs, Cn, Fct, Kmatp, X1v, Fcb); // Form the linear system for particle potential
//     Fcb *= Constants::dt; // Scale the right-hand side vector by the time step

//     Tmatp.reset(Add(1.0, Mmatp, -1.0* Constants::dt, Kmatp));

//     // Solve for the next time-step concentrations
//     Cn.GetTrueDofs(CpV0);
//     Tmatp->Mult(CpV0, RHCp);
//     RHCp += Fcb;

//     Mp_solver.Mult(RHCp, CpVn);

//     // Update only the solid region MAKE INTO FUNCTION
//     for (int p = 0; p < CpV0.Size(); p++){
//         if (PsVc(p) < 1e-5){ // 1e-1 works, but gaps still get smaller
//             (CpVn)(p) = Constants::init_CnC;} // Cp0 initial value

//         if (CpVn(p) < 0.0){
//             (CpVn)(p) = 1.0e-6;} // prevent negative concentrations
//     }

//     Cn.Distribute(CpVn);

//     // Degree of Lithiation
//     utils.CalculateLithiation(Cn, psx, gtPsC);
//     Xfr = utils.GetLithiation();


// }
