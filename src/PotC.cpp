/**
 * @file PotC.cpp
 * @brief Implementation of the potential class for cathode potential simulations.
 */

#include "../include/PotC.hpp"
#include "../include/Constants.hpp"
#include "mfem.hpp"
#include <optional>


PotC::PotC(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc)
    : PotentialBase(geo,para), geometry(geo), domain_parameters(para), boundary_conditions(bc), utils(geo,para), fespace(geo.parfespace), fem(geo.parfespace), dbc_e_bdr(bc.dbc_e_bdr), gtPsi(para.gtPsi), gtPsC(para.gtPsC), 
    ess_tdof_list_e(bc.ess_tdof_list_e), kap(fespace.get()), RpP(fespace.get()), pP0(fespace.get())
    
    {
        cgPP_solver = mfem::CGSolver(MPI_COMM_WORLD);
        B1t = mfem::ParLinearForm(fespace.get());
        X1v = mfem::HypreParVector(fespace.get());
        B1v = mfem::HypreParVector(fespace.get());
        Fpb = mfem::HypreParVector(fespace.get());
        Xs0 = mfem::HypreParVector(fespace.get()); // Initialize the solution vector for particle potential
        RpP = mfem::ParGridFunction(fespace.get());

        Bp2 = std::make_unique<mfem::ParLinearForm>(fespace.get());

        kap = mfem::ParGridFunction(fespace.get()); // Initialize conductivity field
        cKp = mfem::GridFunctionCoefficient(&kap); // Coefficient for conductivity field
        cRp = mfem::GridFunctionCoefficient(&RpP);
        Fpt = mfem::ParLinearForm(fespace.get());

        Kp2 = std::make_unique<mfem::ParBilinearForm>(fespace.get()); // Initialize the bilinear form for conductivity

        pP0 = mfem::ParGridFunction(fespace.get()); // Initialize the potential grid function

        if (gtPsC < 1.0e-200){
            gtPsC = gtPsi;
        }

    }


void PotC::SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{

    BvC = initial_value; // Set the boundary value
    utils.SetInitialValue(ph, BvC); // Initialize potentials

    // kap = psx; // Set the conductivity field to the particle concentration field // HALF
    // kap *= 3.3; // HALF
    fem.InitializeStiffnessMatrix(cKp, Kp2); // Initialize the stiffness matrix
    
    // mfem::ConstantCoefficient dbc_e_Coef(BvC); // Coefficient for Dirichlet boundary conditions // HALF
    // ph.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions // HALF
    
    // fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e); // HALF

    // fem.FormLinearSystem(Kp2, ess_tdof_list_e, ph, B1t, KmP, X1v, B1v); // Assemble the linear system

    // Mpp = std::make_unique<mfem::HypreBoomerAMG>(KmP);  // builds hierarchy once
    // Mpp->SetPrintLevel(0);
    // fem.SolverConditions(KmP, cgPP_solver, *Mpp); // Set up the solver conditions

    // cRp.SetGridFunction(&RpP); // Set the reaction field coefficient
    fem.InitializeForceTerm(cRp, Bp2); // Initialize the force term
    fem.Update(Bp2); // Assemble the force term 
    Fpt = *Bp2; // Assign the force term

    // fem.FormLinearSystem(Kp2, ess_tdof_list_e, ph, Fpt, KmP, X1v, Fpb); // Assemble the force term system HALF

    // mfem::CGSolver cgPP_solver(MPI_COMM_WORLD);
    cgPP_solver.SetRelTol(1e-6);
    cgPP_solver.SetMaxIter(102);
}


void PotC::AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{
	mfem::ConstantCoefficient dbc_e_Coef(BvC);

    ParticleConductivity(Cn, psx); // Update conductivity
    // kap = 1.0; // testing
    fem.Update(Kp2); // Update the stiffness matrix

    potential.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions
    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e); // HALF

    fem.FormLinearSystem(Kp2, ess_tdof_list_e, potential, B1t, KmP, X1v, B1v); // Assemble the linear system

    Mpp = std::make_unique<mfem::HypreBoomerAMG>(KmP);  // builds hierarchy once
    Mpp->SetPrintLevel(0);
    // Mpp->SetOperator(KmP); // Set the operator for the preconditioner
    cgPP_solver.SetPreconditioner(*Mpp);
    cgPP_solver.SetOperator(KmP); // Set the operator for the solver

}

void PotC::AssembleSystemMulti(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups, mfem::ParGridFunction &potential)
{
    mfem::ConstantCoefficient dbc_e_Coef(BvC);

    // Build one effective conductivity field from all particles
    ParticleConductivityMulti(Cn_groups, psi_groups);

    // Update stiffness matrix using the new kap field
    fem.Update(Kp2);

    potential.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr);
    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

    fem.FormLinearSystem(Kp2, ess_tdof_list_e, potential,
                         B1t, KmP, X1v, B1v);

    Mpp = std::make_unique<mfem::HypreBoomerAMG>(KmP);
    Mpp->SetPrintLevel(0);

    cgPP_solver.SetPreconditioner(*Mpp);
    cgPP_solver.SetOperator(KmP);
}


void PotC::UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpP = Rx;
    RpP *= Constants::Frd; // Scale the reaction field

    Bp2->Assemble();
    Fpt = *Bp2; // Assign the force term

    mfem::ConstantCoefficient dbc_e_Coef(BvC);	
    phx.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions
    
    fem.FormLinearSystem(Kp2, ess_tdof_list_e, phx, Fpt, KmP, X1v, Fpb); // Assemble the force term system

    pP0 = phx; // Store the current potential field
    pP0.GetTrueDofs(Xs0); // Extract degrees of freedom

    // std::cout << "PotP Solver" << std::endl;
    cgPP_solver.Mult(Fpb, Xs0); // Solve for the error term
    phx.Distribute(Xs0); // Distribute the updated values

    utils.CalculateGlobalError(pP0, phx, psx, gerror, gtPsC); // Compute global error
    
}



void PotC::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){
    for (int vi = 0; vi < nV; vi++){
        kap(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);
    }
}

void PotC::ParticleConductivityMulti(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups){
    for (int vi = 0; vi < nV; vi++){
        double kap_sum = 0.0;
        for (size_t i = 0; i < Cn_groups.size(); ++i){
            double Cn_val = (*Cn_groups[i])(vi);
            double psx_val = (*psi_groups[i])(vi);
            kap_sum += psx_val * (0.01929 + 0.7045 * tanh(2.399 * Cn_val) - 0.7238 * tanh(2.412 * Cn_val) - 4.2106e-6);
        }
        kap(vi) = kap_sum;
    }
}


