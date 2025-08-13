/**
 * @file PotP.cpp
 * @brief Implementation of the potential class for particle potential simulations.
 */

#include "../include/PotP.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <optional>


// double BvP = 0.0; ///< Global variable for the boundary value of particle potential
// double BvP = -0.1; ///< Global variable for the boundary value of particle potential

PotP::PotP(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_e_bdr(geo.dbc_e_bdr), gtPsi(para.gtPsi), 
    ess_tdof_list_e(geo.ess_tdof_list_e), kap(fespace.get()), RpP(fespace.get()), pP0(fespace.get())
    
    {
    cgPP_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
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
    KmP = std::make_shared<mfem::HypreParMatrix>(); // Initialize the stiffness matrix for conductivity

    pP0 = mfem::ParGridFunction(fespace.get()); // Initialize the potential grid function
    }


void PotP::Initialize(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{


    // dbc_e_bdr = 0; 
    // dbc_e_bdr[0] = 1;
	// mfem::Array<int> ess_tdof_list_e(0);

    BvP = initial_value; // Set the boundary value
    Potentials::SetInitialPotentials(ph, BvP); // Initialize potentials

    kap = psx; // Set the conductivity field to the particle concentration field
    kap *= 3.3;
    // cKp.SetGridFunction(&kap); // Set the conductivity coefficient
    SolverSteps::InitializeStiffnessMatrix(cKp, Kp2); // Initialize the stiffness matrix
    
    mfem::ConstantCoefficient dbc_e_Coef(BvP); // Coefficient for Dirichlet boundary conditions
    ph.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions
    
    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

    SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_e, ph, B1t, KmP, X1v, B1v); // Assemble the linear system
    
    // Mpp.SetOperator(*KmP);
    // mfem::HypreBoomerAMG Mpp(*KmP); // Initialize the preconditioner for the solver
    
    // Mpp.SetPrintLevel(0);
    // Mpp.SetOperator(*KmP);
    // SolverSteps::SolverConditions(KmP, *cgPP_solver, Mpp); // Set up the solver conditions

    Mpp.SetOperator(*KmP);            // build AMG hierarchy once on first matrix
    Mpp.SetPrintLevel(0);

    cgPP_solver->SetRelTol(1e-7);
    cgPP_solver->SetAbsTol(0.0);
    cgPP_solver->SetMaxIter(80);
    cgPP_solver->SetPrintLevel(0);
    // cgPE_solver->SetIterativeMode(true);  // <-- IMPORTANT: use Xe0 as initial guess

    cgPP_solver->SetPreconditioner(Mpp);  // attach once (object stays the same)
    cgPP_solver->SetOperator(*KmP);       // bind initial operator



    // cRp.SetGridFunction(&RpP); // Set the reaction field coefficient
    SolverSteps::InitializeForceTerm(cRp, Bp2); // Initialize the force term
    SolverSteps::Update(Bp2); // Assemble the force term
    // Fpt = std::move(*Bp2); // Move the force term
    Fpt = *Bp2; // Assign the force term

    SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_e, ph, Fpt, KmP, X1v, Fpb); // Assemble the force term system
}


void PotP::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{


    // // dbc_e_bdr = 0; 
    // // dbc_e_bdr[0] = 1;
	// // mfem::Array<int> ess_tdof_list_e(0);

    // // ParticleConductivity(Cn, psx); // Update conductivity
    // // cKp.SetGridFunction(&kap); // Set the conductivity coefficient

    // // SolverSteps::Update(Kp2); // Update the stiffness matrix
    // mfem::ConstantCoefficient dbc_e_Coef(BvP); // Coefficient for Dirichlet boundary conditions
    // potential.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions
    
    // fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);
    // // std::cout << "[PotP::TimeStep] Essential DOFs: " << ess_tdof_list_e.Size() << std::endl;
    
    // // SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_e, potential, B1t, KmP, X1v, B1v); // Assemble the linear system

    // // Mpp.SetOperator(*KmP); // Set the preconditioner operator
    // // cgPP_solver->SetPreconditioner(Mpp); // Attach the preconditioner to the solver
    // // cgPP_solver->SetOperator(*KmP); // Set the operator for the solver 

	mfem::ConstantCoefficient dbc_e_Coef(BvP);	



}

void PotP::Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpP = Rx;
    RpP *= Constants::Frd; // Scale the reaction field

    Bp2 -> Assemble();
    Fpt = *Bp2; // Assign the force term

    cout << "BvP: " << BvP << endl;
    mfem::ConstantCoefficient dbc_e_Coef(BvP);	

    phx.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions
    // Kp2->FormLinearSystem(ess_tdof_list_e, phx, Fpt, *KmP, X1v, Fpb); // Assemble the force term system

    SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_e, phx, Fpt, KmP, X1v, Fpb); // Assemble the force term system


    pP0 = phx; // Store the current potential field
    pP0.GetTrueDofs(Xs0); // Extract degrees of freedom

    // Rebind operator & preconditioner to this KmP **every time**
    Mpp.SetOperator(*KmP);
    // cgPP_solver->SetPreconditioner(Mpp);
    cgPP_solver->SetOperator(*KmP);
    // cgPP_solver->iterative_mode = false;

    cgPP_solver->Mult(Fpb, Xs0); // Solve for the error term
    phx.Distribute(Xs0); // Distribute the updated values

    Potentials::ComputeGlobalError(pP0, phx, psx, gerror, gtPsi); // Compute global error
    
    
    
    
    
    
    
    
    // dbc_e_bdr = 0; 
    // dbc_e_bdr[0] = 1;
	// mfem::Array<int> ess_tdof_list_e(0);

    // Potentials::AssembleForceVector(Rx, RpP, Constants::Frd, cRp, Bp2, Fpt); // Create reaction field
    
    // mfem::ConstantCoefficient dbc_e_Coef(BvP); // Coefficient for Dirichlet boundary conditions
    // phx.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions
    
    // fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);
    // // std::cout << "[PotP::Advance] Essential DOFs: " << ess_tdof_list_e.Size() << std::endl;
    
    // SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_e, phx, Fpt, KmP, X1v, Fpb); // Assemble the force term system
    // // Kp2->FormLinearSystem(ess_tdof_list_e, phx, Fpt, *KmP, X1v, Fpb);


    // pP0 = phx; // Store the current potential field
    // pP0.GetTrueDofs(Xs0); // Extract degrees of freedom
    
    // Mpp.SetOperator(*KmP); // Set the preconditioner operator
    // // Mpp.SetPrintLevel(0);

    // // cgPP_solver->SetPreconditioner(Mpp); // Attach the preconditioner to the solver
    // // cgPP_solver->SetOperator(*KmP); // Set the operator for the solver 
    
    // SolverSteps::SolverConditions(KmP, *cgPP_solver, Mpp); // Set up the solver conditions

    // cgPP_solver->Mult(Fpb, Xs0); // Solve for the error term
    // phx.Distribute(Xs0); // Distribute the updated values

    // Potentials::ComputeGlobalError(pP0, phx, psx, gerror, gtPsi); // Compute global error
}






// void PotP::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
// {
//     ParticleConductivity(Cn, psx); // Update conductivity
//     auto cKp = std::make_shared<mfem::GridFunctionCoefficient>(kap.get());

//     Potentials::ImplementBoundaryConditions(dbc_e_Coef, BvP, potential, dbc_e_bdr); // Apply BCs
//     SolverSteps::StiffnessMatrix(cKp, ess_tdof_list_e, potential, B1t, KmP, B1v);
//     Potentials::PCG_Solver(Mpp, *cgPP_solver, *KmP); // Solve the system
// }

void PotP::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){
    for (int vi = 0; vi < nV; vi++){
        kap(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);
        // kap(vi) = 3.3 * psx(vi);
    }
}

// void PotP::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) 
// {
//     Potentials::CreateReaction(Rx, *RpP, Constants::Frd); // Create reaction field

//     // Assemble the force term without applying boundary conditions
//     SolverSteps::ForceTerm(*RpP, ftPotP); // false since not applying BCs

//     Potentials::ForceVector(*K, ess_tdof_list_e, phx, ftPotP, *KmP, X1v, Fpb, dbc_e_Coef, dbc_e_bdr); // Force vector
//     // Potentials::ErrorCalculation(phx, *cgPP_solver, Fpb, psx, error_P, gerror, gtPsi); // Compute global error
// }


