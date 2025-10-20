/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "../include/PotE.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include "../include/CnE.hpp"
#include <optional>

PotE::PotE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_bdr(geo.dbc_bdr), gtPse(para.gtPse), 
    kpl(fespace.get()), RpE(fespace.get()), Dmp(fespace.get()), pE0(fespace.get()), ess_tdof_potE(geo.ess_tdof_potE), phE_bc(fespace.get())
    
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

#include <mpi.h>

void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{


    const int anchor = 2000;

    if (!anchor_set) {
        ess_tdof_potE.SetSize(1);
        ess_tdof_potE[0] = anchor;  
        anchor_set = true; 
    }


    
    BvE = initial_value; // Set the boundary value.
    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials

    // std::cout << "PotE: Ess tdof for potential on electrode: " << ess_tdof_potE[0] << std::endl;

    phE_bc = ph;                       // start from current iterate
    if (ess_tdof_potE.Size() == 1) {
        mfem::Vector td; 
        phE_bc.GetTrueDofs(td); // get true dofs
        td(ess_tdof_potE[0]) = Constants::init_BvE;     // set anchor value
        phE_bc.SetFromTrueDofs(td); // set grid function from modified true dofs
        ph = phE_bc;
    }
    
    SolverSteps::InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix

    // mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    // ph.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF

    // fespace->GetEssentialTrueDofs(dbc_bdr, ess_tdof_list_potE); // Get essential true degrees of freedom for Dirichlet boundary conditions HALF

    // std::cout << "Number of essential true dofs for PotE: " << ess_tdof_list_potE.Size() << std::endl;
    // std::cout << "boundary_dofs size for PotE: " << boundary_dofs.Size() << std::endl;

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system HALF
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system FULL

    std::cout << "Linear system for PotE assembled." << std::endl;

    Mpe = std::make_unique<mfem::HypreBoomerAMG>(Kml);  // builds hierarchy once
    Mpe->SetPrintLevel(0);
    SolverSteps::SolverConditions(Kml, cgPE_solver, *Mpe); // Set up the solver conditions

    std::cout << "Solver conditions for PotE set." << std::endl;

    SolverSteps::InitializeForceTerm(cRe, Bl2); // Initialize the force term
    SolverSteps::Update(Bl2); // Update the force term
    Flt = *Bl2; // Move the force term

    std::cout << "Force term for PotE initialized." << std::endl;

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, Flt, Kml, X1v, Flb); // Assemble the force term system HALF
    SolverSteps::FormLinearSystem(Kl2, boundary_dofs, ph, Flt, Kml, X1v, Flb); // Assemble the force term system FULL

    std::cout << "Force term system for PotE assembled." << std::endl;

    SolverSteps::InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, ph, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    std::cout << "Diffusivity matrix for PotE initialized." << std::endl;
 

}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn)
{
    ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    SolverSteps::Update(Kl1); // Update the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm.Mult(CeVn, LpCe); 

    SolverSteps::Update(Kl2); // Update the conductivity matrix

    // mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    // potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system HALF
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system FULL

    Mpe->SetOperator(Kml); // Set the operator for the preconditioner
    cgPE_solver.SetPreconditioner(*Mpe);
    cgPE_solver.SetOperator(Kml); // Set the operator for the solver

}

void PotE::Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpE = Rx;
    RpE.Neg();

    Bl2->Assemble();
    Flt = *Bl2;

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    phx.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions
    
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
    RpE = Rx1;
    RpE += Rx2;

    RpE.Neg(); // check and see if this should be negative
    
    Bl2->Assemble();
    Flt = *Bl2;

    // // Enforce the same gauge
    // phE_bc = phx;
    // if (ess_tdof_potE.Size() == 1) {
    //     mfem::Vector td; phE_bc.GetTrueDofs(td);
    //     td(ess_tdof_potE[0]) = BvE;
    //     phE_bc.SetFromTrueDofs(td);
    //     phx = phE_bc;
    // }

    phx(ess_tdof_potE[0]) = BvE; // Enforce the gauge directly on the solution vector
    // std::cout << "Number of essential true dofs for PotE: " << ess_tdof_potE.Size() << std::endl;
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, phx, Flt, Kml, X1v, Flb);

    RHSl = Flb;
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values
    Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
}




void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }
}
