/**
 * @file PotA.cpp
 * @brief Implementation of the potential class for anode potential simulations.
 */

#include "../include/PotA.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <optional>

PotA::PotA(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_w_bdr(geo.dbc_w_bdr), gtPsi(para.gtPsi), 
    ess_tdof_list_w(geo.ess_tdof_list_w), kap(fespace.get()), RpP(fespace.get()), pP0(fespace.get()), gtPsA(para.gtPsA)
    
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

    if (gtPsA < 1.0e-200){
        gtPsA = gtPsi;
    }

    }


void PotA::Initialize(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{

    BvA = initial_value; // Set the boundary value
    Potentials::SetInitialPotentials(ph, BvA); // Initialize potentials

    kap = psx; // Set the conductivity field to the particle concentration field
    kap *= 3.3;
    SolverSteps::InitializeStiffnessMatrix(cKp, Kp2); // Initialize the stiffness matrix
    
    mfem::ConstantCoefficient dbc_w_Coef(BvA); // Coefficient for Dirichlet boundary conditions
    ph.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); // Apply Dirichlet boundary conditions
    
    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

    SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_w, ph, B1t, KmP, X1v, B1v); // Assemble the linear system

    Mpp = std::make_unique<mfem::HypreBoomerAMG>(KmP);  // builds hierarchy once
    Mpp->SetPrintLevel(0);
    SolverSteps::SolverConditions(KmP, cgPP_solver, *Mpp); // Set up the solver conditions

    // cRp.SetGridFunction(&RpP); // Set the reaction field coefficient
    SolverSteps::InitializeForceTerm(cRp, Bp2); // Initialize the force term
    SolverSteps::Update(Bp2); // Assemble the force term
    Fpt = *Bp2; // Assign the force term

    SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_w, ph, Fpt, KmP, X1v, Fpb); // Assemble the force term system
}


void PotA::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{
	mfem::ConstantCoefficient dbc_w_Coef(BvA);

}

void PotA::Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpP = Rx;
    RpP *= Constants::Frd; // Scale the reaction field

    Bp2->Assemble();
    Fpt = *Bp2; // Assign the force term

    mfem::ConstantCoefficient dbc_w_Coef(BvA);	

    phx.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); // Apply Dirichlet boundary conditions
    
    SolverSteps::FormLinearSystem(Kp2, ess_tdof_list_w, phx, Fpt, KmP, X1v, Fpb); // Assemble the force term system

    pP0 = phx; // Store the current potential field
    pP0.GetTrueDofs(Xs0); // Extract degrees of freedom

    cgPP_solver.Mult(Fpb, Xs0); // Solve for the error term
    phx.Distribute(Xs0); // Distribute the updated values

    Potentials::ComputeGlobalError(pP0, phx, psx, gerror, gtPsi); // Compute global error
    
}



void PotA::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){
    for (int vi = 0; vi < nV; vi++){
        kap(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);
    }
}



