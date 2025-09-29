/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "../include/PotE.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include "../include/CnE.hpp"
#include <optional>


// double BvE = 0.0; ///< Global variable for the boundary value of electrolyte potential
// double BvE = -0.4686; ///< Global variable for the boundary value of electrolyte potential

PotE::PotE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_CnE_bdr(geo.dbc_CnE_bdr), gtPse(para.gtPse), 
    kpl(fespace.get()), RpE(fespace.get()), Dmp(fespace.get()), pE0(fespace.get())
    
    {

    cgPE_solver = mfem::CGSolver(MPI_COMM_WORLD);
    
    B1t = mfem::ParLinearForm(fespace.get());
    X1v = mfem::HypreParVector(fespace.get());
    B1v = mfem::HypreParVector(fespace.get());
    Flb = mfem::HypreParVector(fespace.get());
    LpCe = mfem::HypreParVector(fespace.get()); // Initialize the vector for concentration degrees of freedom
    RpE = mfem::ParGridFunction(fespace.get());
    Dmp = mfem::ParGridFunction(fespace.get()); // Initialize diffusivity field

    Bl2 = std::make_unique<mfem::ParLinearForm>(fespace.get());

    kpl = mfem::ParGridFunction(fespace.get()); // Initialize conductivity field
    cKe = mfem::GridFunctionCoefficient(&kpl); // Coefficient for conductivity field
    cRe = mfem::GridFunctionCoefficient(&RpE);
    Flt = mfem::ParLinearForm(fespace.get());
    cDm = mfem::GridFunctionCoefficient(&Dmp); // Coefficient for diffusivity field

    Kl1 = std::make_unique<mfem::ParBilinearForm>(fespace.get()); // Initialize the bilinear form for conductivity
    Kl2 = std::make_unique<mfem::ParBilinearForm>(fespace.get()); // Initialize the bilinear form for conductivity

    pE0 = mfem::ParGridFunction(fespace.get()); // Initialize the potential grid function
    Xe0 = mfem::HypreParVector(fespace.get()); // Initialize the solution vector for potential
    RHSl = mfem::HypreParVector(fespace.get()); // Initialize the right-hand side vector for potential

    }

void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{
    BvE = initial_value; // Set the boundary value.
    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials

    SolverSteps::InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    ph.ProjectBdrCoefficient(dbc_potE_Coef, dbc_CnE_bdr); // Apply Dirichlet boundary conditions 

    fespace->GetEssentialTrueDofs(dbc_CnE_bdr, ess_tdof_list_potE); // Get essential true degrees of freedom for Dirichlet boundary conditions

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system

    Mpe = std::make_unique<mfem::HypreBoomerAMG>(Kml);  // builds hierarchy once
    Mpe->SetPrintLevel(0);
    SolverSteps::SolverConditions(Kml, cgPE_solver, *Mpe); // Set up the solver conditions

    SolverSteps::InitializeForceTerm(cRe, Bl2); // Initialize the force term
    SolverSteps::Update(Bl2); // Update the force term
    Flt = *Bl2; // Move the force term

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, Flt, Kml, X1v, Flb); // Assemble the force term system

    SolverSteps::InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, ph, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system
 

}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn)
{
    ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    SolverSteps::Update(Kl1); // Update the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm.Mult(CeVn, LpCe); 

    SolverSteps::Update(Kl2); // Update the conductivity matrix

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_CnE_bdr); // Apply Dirichlet boundary conditions

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system

    Mpe->SetOperator(Kml); // Set the operator for the preconditioner
    cgPE_solver.SetPreconditioner(*Mpe);
    cgPE_solver.SetOperator(Kml); // Set the operator for the solver

}

void PotE::Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    // Potentials::AssembleForceVector(Rx, RpE, -1.0, cRe, Bl2, Flt); // Create reaction field
    RpE = Rx;
    RpE.Neg();

    Bl2->Assemble();
    Flt = *Bl2;

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    phx.ProjectBdrCoefficient(dbc_potE_Coef, dbc_CnE_bdr); // Apply Dirichlet boundary conditions
    
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, phx, Flt, Kml, X1v, Flb); // Assemble the force term system

    RHSl = Flb;
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values

    Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
}

void PotE::Advance(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    // Potentials::AssembleForceVector(Rx, RpE, -1.0, cRe, Bl2, Flt); // Create reaction field
    Rx2 += Rx1;
    RpE = Rx2;
    RpE.Neg();

    Bl2->Assemble();
    Flt = *Bl2;

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    phx.ProjectBdrCoefficient(dbc_potE_Coef, dbc_CnE_bdr); // Apply Dirichlet boundary conditions
    
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, phx, Flt, Kml, X1v, Flb); // Assemble the force term system

    RHSl = Flb;
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values

    // Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
}


void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }
}
