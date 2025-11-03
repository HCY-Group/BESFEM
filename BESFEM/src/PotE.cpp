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
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_bdr(geo.dbc_bdr), gtPse(para.gtPse), dbc_e_bdr(geo.ess_tdof_marker),
    kpl(fespace.get()), RpE(fespace.get()), Dmp(fespace.get()), pE0(fespace.get()), ess_tdof_potE(geo.ess_tdof_listPinned), phE_bc(fespace.get()), CeVn(fespace.get())
    
    {
    cgPE_solver = mfem::CGSolver(MPI_COMM_WORLD); // Initialize the conjugate gradient solver for potential
    
    B1t = mfem::ParLinearForm(fespace.get()); // Initialize the linear form for the potential
    X1v = mfem::HypreParVector(fespace.get()); // Initialize the solution vector for potential
    B1v = mfem::HypreParVector(fespace.get()); // Initialize the right-hand side vector for potential
    Flb = mfem::HypreParVector(fespace.get()); // Initialize the force term vector for potential
    LpCe = mfem::HypreParVector(fespace.get()); // Initialize the vector for concentration degrees of freedom
    RpE = mfem::ParGridFunction(fespace.get()); // Initialize reaction rate field
    Dmp = mfem::ParGridFunction(fespace.get()); // Initialize diffusivity field

    // Bl2 = std::make_unique<mfem::ParLinearForm>(fespace.get());

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
    BvE = initial_value; // Set the boundary value.
    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials

    SolverSteps::InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix

    // if (mfem::Mpi::WorldRank() == geometry.rkpp && ess_tdof_potE.Size() > 0)
    // {
    //     ph(ess_tdof_potE[0]) = BvE;
    // }

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system FULL

    // std::cout << " finish form linear system PotE " << std::endl;

    Mpe = std::make_unique<mfem::HypreBoomerAMG>(Kml);  // builds hierarchy once
    Mpe->SetPrintLevel(0);
    cgPE_solver.SetRelTol(1e-6);
    cgPE_solver.SetMaxIter(102);
    cgPE_solver.SetPreconditioner(*Mpe);
    cgPE_solver.SetOperator(Kml); 

    // SolverSteps::SolverConditions(Kml, cgPE_solver, *Mpe); // Set up the solver conditions


    // Mpp = std::make_unique<mfem::HypreBoomerAMG>(KmP);  // builds hierarchy once
    // Mpp->SetPrintLevel(0);
    // cgPP_solver.SetRelTol(1e-6);
    // cgPP_solver.SetMaxIter(102);
    // cgPP_solver.SetPreconditioner(*Mpp);
    // cgPP_solver.SetOperator(KmP); 

    // SolverSteps::SolverConditions(KmP, cgPP_solver, *Mpp); // Set up the solver conditions

    // std::cout << " finish solver conditions PotE " << std::endl;

    Bl2 = std::make_unique<mfem::ParLinearForm>(fespace.get());
    Bl2->AddDomainIntegrator(new mfem::DomainLFIntegrator(cRe));
    // SolverSteps::InitializeForceTerm(cRe, Bl2); // Initialize the force term
    // std::unique_ptr<mfem::ParLinearForm>Bl2(new mfem::ParLinearForm(fespace.get()));
    // Bl2->AddDomainIntegrator(new mfem::DomainLFIntegrator(cRe));
    // Bl2->Update();
    Bl2->Assemble();
    // SolverSteps::Update(Bl2); // Update the force term
    Flt = *Bl2; // Move the force term
    // Flt = std::move(*Bl2);

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, Flt, Kml, X1v, Flb); // Assemble the force term system HALF

    SolverSteps::FormLinearSystem(Kl2, boundary_dofs, ph, Flt, Kml, X1v, Flb); // Assemble the force term system FULL

    SolverSteps::InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, ph, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system
 
    // std::cout << " exit Initialize PotE " << std::endl;

}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{

    // mfem::ConstantCoefficient dbc_e_Coef(BvE);

    // ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }

    // SolverSteps::Update(Kl1); // Update the diffusivity matrix
    cDm.SetGridFunction(&Dmp);
    Kl1->Update();
    Kl1->Assemble();
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    // std::cout << "boundary_dofs size: " << boundary_dofs.Size() << " on rank " << myid << std::endl;

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm.Mult(CeVn, LpCe); 

    // mfem::HypreParVector LpCe_true(fespace.get());
    // Kdm.Mult(CeVn, LpCe_true);
    // LpCe = LpCe_true;


    // SolverSteps::Update(Kl2); // Update the conductivity matrix
    cKe.SetGridFunction(&kpl);
    Kl2->Update();
    Kl2->Assemble();

    // if (mfem::Mpi::WorldRank() == geometry.rkpp) {
    //     // std::cout << "pinning on rank: " << mfem::Mpi::WorldRank() << std::endl;
    //     potential(ess_tdof_potE[0]) = BvE;
    // }

    // mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    // potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF
    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system HALF

    // mfem::ConstantCoefficient dbc_e_Coef(BvE);


    // potential.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); // Apply Dirichlet boundary conditions

    // if (mfem::Mpi::WorldRank() == geometry.rkpp) {
    //     // std::cout << "pinning on rank: " << mfem::Mpi::WorldRank() << std::endl;
    //     potential(ess_tdof_potE[0]) = BvE;
    //     // phx.Save("phE_pin_2");
    // } 
    
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, potential, Flt, Kml, X1v, Flb); // Assemble the conductivity matrix system FULL


    Mpe->SetOperator(Kml); // Set the operator for the preconditioner
    cgPE_solver.SetPreconditioner(*Mpe);
    cgPE_solver.SetOperator(Kml); // Set the operator for the solver

    // std::cout << " exit TimeStep PotE " << std::endl;

}


void PotE::Advance(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{    
    RpE = Rx2;
    RpE += Rx1;
    RpE.Neg(); // check and see if this should be negative
    
    // SolverSteps::InitializeForceTerm(cRe, Bl2); // Initialize the force term
    Bl2 = std::make_unique<mfem::ParLinearForm>(fespace.get());
    Bl2->AddDomainIntegrator(new mfem::DomainLFIntegrator(cRe));
    // Bl2->Update();
    Bl2->Assemble();
    // SolverSteps::Update(Bl2); // Update the force term
    Flt = *Bl2; // Move the force term


    // if (mfem::Mpi::WorldRank() == geometry.rkpp) {
    //     phx(ess_tdof_potE[0]) = BvE;
    // } 

    
    Kl2->Update();
    Kl2->Assemble();

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, phx, Flt, Kml, X1v, Flb);
    // Kl2->FormLinearSystem(ess_tdof_potE, phx, Flt, Kml, X1v, Flb);
    Kl2->FormLinearSystem(ess_tdof_potE, phx, Flt, Kml, X1v, Flb);

    RHSl = 0.0;
    RHSl = Flb; // when this is commented out, same results in serial and parallel
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom
    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values
    Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error

    // double ref_val = phx.Sum() / phx.Size();
    // phx -= ref_val;


}



























void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn)
{
    myid = mfem::Mpi::WorldRank();


    // ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }

    // SolverSteps::Update(Kl1); // Update the diffusivity matrix
    cDm.SetGridFunction(&Dmp);
    Kl1->Update();
    Kl1->Assemble();
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    // std::cout << "boundary_dofs size: " << boundary_dofs.Size() << " on rank " << myid << std::endl;

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm.Mult(CeVn, LpCe); 

    // mfem::HypreParVector LpCe_true(fespace.get());
    // Kdm.Mult(CeVn, LpCe_true);
    // LpCe = LpCe_true;


    std::cout << "[Rank " << myid << "] LpCe norm before BCs: " << LpCe.Norml2() << std::endl;



    // SolverSteps::Update(Kl2); // Update the conductivity matrix
    cKe.SetGridFunction(&kpl);
    Kl2->Update();
    Kl2->Assemble();

    // if (mfem::Mpi::WorldRank() == geometry.rkpp) {
    //     // std::cout << "pinning on rank: " << mfem::Mpi::WorldRank() << std::endl;
    //     potential(ess_tdof_potE[0]) = BvE;
    // }

    // mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    // potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF
    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system HALF

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, potential, Flt, Kml, X1v, Flb); // Assemble the conductivity matrix system FULL


    Mpe->SetOperator(Kml); // Set the operator for the preconditioner
    cgPE_solver.SetPreconditioner(*Mpe);
    cgPE_solver.SetOperator(Kml); // Set the operator for the solver

    // std::cout << " exit TimeStep PotE " << std::endl;

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

void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }
}
