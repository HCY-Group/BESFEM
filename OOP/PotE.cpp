/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "PotE.hpp"
#include "../code/Constants.hpp"
#include "mfem.hpp"
#include "CnE.hpp"
#include <optional>


// double BvE = 0.0; ///< Global variable for the boundary value of electrolyte potential
// double BvE = -0.4686; ///< Global variable for the boundary value of electrolyte potential

PotE::PotE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_w_bdr(geo.dbc_w_bdr), gtPse(para.gtPse), 
    ess_tdof_list_w(geo.ess_tdof_list_w), kpl(fespace.get()), RpE(fespace.get()), Dmp(fespace.get()), pE0(fespace.get())
    
    {

    cgPE_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
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
    Kml = std::make_shared<mfem::HypreParMatrix>(); // Initialize the stiffness matrix for conductivity
    Kdm = std::make_shared<mfem::HypreParMatrix>(); // Initialize the stiffness matrix for conductivity

    pE0 = mfem::ParGridFunction(fespace.get()); // Initialize the potential grid function
    Xe0 = mfem::HypreParVector(fespace.get()); // Initialize the solution vector for potential
    RHSl = mfem::HypreParVector(fespace.get()); // Initialize the right-hand side vector for potential

    }

void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{
    BvE = initial_value; // Set the boundary value.
    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials

    SolverSteps::InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix

    mfem::ConstantCoefficient dbc_w_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    ph.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); // Apply Dirichlet boundary conditions 

    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w); // Get essential true degrees of freedom for Dirichlet boundary conditions
    // std::cout << "[PotE::Initialize] Essential DOFs: " << ess_tdof_list_w.Size() << std::endl;

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_w, ph, B1t, Kml, X1v, B1v); // Assemble the linear system
    
    mfem::HypreBoomerAMG Mpe(*Kml);
    SolverSteps::SolverConditions(Kml, *cgPE_solver, Mpe); // Set up the solver conditions

    SolverSteps::InitializeForceTerm(cRe, Bl2); // Initialize the force term
    Flt = std::move(*Bl2); // Move the force term

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_w, ph, Flt, Kml, X1v, Flb); // Assemble the force term system

    SolverSteps::InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, ph, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn)
{
    ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity
    cDm.SetGridFunction(&Dmp); // Set the diffusivity coefficient

    SolverSteps::Update(Kl1); // Update the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm->Mult(CeVn, LpCe); // Multiply concentration by diffus

    cKe.SetGridFunction(&kpl); // Set the conductivity coefficient
    SolverSteps::Update(Kl2); // Update the conductivity matrix

    mfem::ConstantCoefficient dbc_w_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    potential.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); // Apply Dirichlet boundary conditions

    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w); // Get essential true degrees of freedom for Dirichlet boundary conditions

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_w, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system

    Mpe.SetOperator(*Kml); // Set the operator for the preconditioner
    cgPE_solver->SetPreconditioner(Mpe); // Attach the preconditioner to the solver
    cgPE_solver->SetOperator(*Kml); // Set the operator for the solver

}

void PotE::Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{

    Potentials::AssembleForceVector(Rx, RpE, -1.0, cRe, Bl2, Flt); // Create reaction field
    
    mfem::ConstantCoefficient dbc_w_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    phx.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); // Apply Dirichlet boundary conditions
    
    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w); // Get essential true degrees of freedom for Dirichlet boundary conditions
    
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_w, phx, Flt, Kml, X1v, Flb); // Assemble the force term system

    RHSl = Flb;
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom
    Mpe.SetOperator(*Kml); // Set the preconditioner operator
    cgPE_solver->SetPreconditioner(Mpe); // Attach the preconditioner to the solver
    cgPE_solver->SetOperator(*Kml); // Set the operator for the solver 
    cgPE_solver->Mult(RHSl, Xe0); // Solve for the error term
    phx.Distribute(Xe0); // Distribute the updated values

    Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
}

// void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
// {
//     ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity
//     auto cDm = std::make_shared<mfem::GridFunctionCoefficient>(Dmp.get());
//     auto cKe = std::make_shared<mfem::GridFunctionCoefficient>(kpl.get());

//     SolverSteps::StiffnessMatrix(cDm, boundary_dofs, potential, B1t, Kdm, B1v);
//     Cn.GetTrueDofs(*CeVn); 
//     Kdm->Mult(*CeVn, *LpCe); // Multiply concentration by diffusivity matrix
    
//     Potentials::ImplementBoundaryConditions(dbc_w_Coef, BvE, potential, dbc_w_bdr); // Apply boundary conditions
//     SolverSteps::StiffnessMatrix(cKe, ess_tdof_list_w, potential, B1t, KmE, B1v);
//     Potentials::PCG_Solver(Mpe, *cgPE_solver, *KmE); // Solve the system
// }

void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }
}

// void PotE::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) 
// {
//     Potentials::CreateReaction(Rx, *RpE, -1.0); // Create reaction field

//     // Assemble the force term without applying boundary conditions
//     SolverSteps::ForceTerm(*RpE, ftPotE); // false since not applying BCs
    
//     Potentials::ForceVector(*K, ess_tdof_list_w, phx, ftPotE, *KmE, X1v, Flb, dbc_w_Coef, dbc_w_bdr); // Force vector
//     RHSl = Flb;
//     RHSl += *LpCe;
//     Potentials::ErrorCalculation(phx, *cgPE_solver, RHSl, psx, error_E, gerror, gtPse); // Compute global error
// }