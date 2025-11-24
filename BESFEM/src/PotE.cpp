/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "../include/PotE.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include "../include/CnE.hpp"
#include <optional>
#include "../include/SimTypes.hpp"
#include <mpi.h>


// using sim::CellMode;

PotE::PotE(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::CellMode mode)
    : PotentialBase(geo,para), geometry(geo), domain_parameters(para), boundary_conditions(bc), utils(geo,para), fespace(geo.parfespace), fem(geo.parfespace), dbc_bdr(bc.dbc_bdr), gtPse(para.gtPse),
    kpl(fespace.get()), RpE(fespace.get()), Dmp(fespace.get()), pE0(fespace.get()), ess_tdof_potE(bc.ess_tdof_list), phE_bc(fespace.get()), CeVn(fespace.get()), mode_(mode)
    
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


void PotE::SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{
    if(mode_ == sim::CellMode::HALF){
        std::cout << "Initializing PotE for Half Cell" << std::endl;

        BvE = initial_value; // Set the boundary value.
        utils.SetInitialValue(ph, BvE); // Initialize potentials

        fem.InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix

        mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
        ph.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions 

        fespace->GetEssentialTrueDofs(dbc_bdr, ess_tdof_potE); // Get essential true degrees of freedom for Dirichlet boundary conditions

        fem.FormLinearSystem(Kl2, ess_tdof_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system

        Mpe = std::make_unique<mfem::HypreBoomerAMG>(Kml);  // builds hierarchy once
        Mpe->SetPrintLevel(0);
        fem.SolverConditions(Kml, cgPE_solver, *Mpe); // Set up the solver conditions

        fem.InitializeForceTerm(cRe, Bl2); // Initialize the force term
        fem.Update(Bl2); // Update the force term
        Flt = *Bl2; // Move the force term

        fem.FormLinearSystem(Kl2, ess_tdof_potE, ph, Flt, Kml, X1v, Flb); // Assemble the force term system

        fem.InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
        fem.FormLinearSystem(Kl1, boundary_dofs, ph, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system
    }

    else if (mode_ == sim::CellMode::FULL){
    
        BvE = initial_value; // Set the boundary value.
        utils.SetInitialValue(ph, BvE); // Initialize potentials

        fem.InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix
        fem.FormLinearSystem(Kl2, ess_tdof_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system FULL

        Mpe = std::make_unique<mfem::HypreBoomerAMG>(Kml);  // builds hierarchy once
        Mpe->SetPrintLevel(0);
        cgPE_solver.SetRelTol(1e-6);
        cgPE_solver.SetMaxIter(102);
        cgPE_solver.SetPreconditioner(*Mpe);
        cgPE_solver.SetOperator(Kml); 

        Bl2 = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Bl2->AddDomainIntegrator(new mfem::DomainLFIntegrator(cRe));
        Bl2->Assemble();
        Flt = *Bl2; // Move the force term

        fem.FormLinearSystem(Kl2, boundary_dofs, ph, Flt, Kml, X1v, Flb); // Assemble the force term system FULL

        fem.InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
        fem.FormLinearSystem(Kl1, boundary_dofs, ph, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system
    }
}

void PotE::AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{
    if(mode_ == sim::CellMode::HALF){

        for (int vi = 0; vi < nV; vi++){
            dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
            Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
            kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
        }

        fem.Update(Kl1); // Update the diffusivity matrix
        fem.FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

        Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
        Kdm.Mult(CeVn, LpCe); 

        fem.Update(Kl2); // Update the conductivity matrix

        mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
        potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions

        fem.FormLinearSystem(Kl2, ess_tdof_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system

        Mpe->SetOperator(Kml); // Set the operator for the preconditioner
        cgPE_solver.SetPreconditioner(*Mpe);
        cgPE_solver.SetOperator(Kml); // Set the operator for the solver

    }

    if(mode_ == sim::CellMode::FULL){

        for (int vi = 0; vi < nV; vi++){
            dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
            Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
            kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
        }

        cDm.SetGridFunction(&Dmp);
        Kl1->Update();
        Kl1->Assemble();
        fem.FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

        Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
        Kdm.Mult(CeVn, LpCe); 

        cKe.SetGridFunction(&kpl);
        Kl2->Update();
        Kl2->Assemble();

        fem.FormLinearSystem(Kl2, ess_tdof_potE, potential, Flt, Kml, X1v, Flb); // Assemble the conductivity matrix system FULL

        Mpe->SetOperator(Kml); // Set the operator for the preconditioner
        cgPE_solver.SetPreconditioner(*Mpe);
        cgPE_solver.SetOperator(Kml); // Set the operator for the solver
    }

}


void PotE::UpdatePotential(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{  
    RpE = Rx2;
    RpE += Rx1;
    RpE.Neg(); // check and see if this should be negative
    
    Bl2 = std::make_unique<mfem::ParLinearForm>(fespace.get());
    Bl2->AddDomainIntegrator(new mfem::DomainLFIntegrator(cRe));
    Bl2->Assemble();
    Flt = *Bl2; // Move the force term

    
    Kl2->Update();
    Kl2->Assemble();

    Kl2->FormLinearSystem(ess_tdof_potE, phx, Flt, Kml, X1v, Flb);

    RHSl = 0.0;
    RHSl = Flb; // when this is commented out, same results in serial and parallel
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom
    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values
    utils.CalculateGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error

}

void PotE::UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpE = Rx;
    RpE.Neg();

    Bl2->Assemble();
    Flt = *Bl2;

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    phx.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions
    
    fem.FormLinearSystem(Kl2, ess_tdof_potE, phx, Flt, Kml, X1v, Flb); // Assemble the force term system

    RHSl = Flb;
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values

    utils.CalculateGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
}

void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }
}

































void PotE::AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn)
{
    myid = mfem::Mpi::WorldRank();


    // ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }

    // fem.Update(Kl1); // Update the diffusivity matrix
    cDm.SetGridFunction(&Dmp);
    Kl1->Update();
    Kl1->Assemble();
    fem.FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    // std::cout << "boundary_dofs size: " << boundary_dofs.Size() << " on rank " << myid << std::endl;

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm.Mult(CeVn, LpCe); 


    // fem.Update(Kl2); // Update the conductivity matrix
    cKe.SetGridFunction(&kpl);
    Kl2->Update();
    Kl2->Assemble();

    // if (mfem::Mpi::WorldRank() == geometry.rkpp) {
    //     // std::cout << "pinning on rank: " << mfem::Mpi::WorldRank() << std::endl;
    //     potential(ess_tdof_potE[0]) = BvE;
    // }

    // mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    // potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF
    // fem.FormLinearSystem(Kl2, ess_tdof_list_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system HALF

    fem.FormLinearSystem(Kl2, ess_tdof_potE, potential, Flt, Kml, X1v, Flb); // Assemble the conductivity matrix system FULL


    Mpe->SetOperator(Kml); // Set the operator for the preconditioner
    cgPE_solver.SetPreconditioner(*Mpe);
    cgPE_solver.SetOperator(Kml); // Set the operator for the solver

}


// void PotE::Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
// {
//     RpE = Rx;
//     RpE.Neg();

//     Bl2->Assemble();
//     Flt = *Bl2;

//     mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
//     phx.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions
    
//     fem.FormLinearSystem(Kl2, ess_tdof_list_potE, phx, Flt, Kml, X1v, Flb); // Assemble the force term system


//     RHSl = Flb;
//     RHSl += LpCe;

//     pE0 = phx; // Store the current potential field
//     pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

//     cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

//     phx.Distribute(Xe0); // Distribute the updated values

//     Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
// }

// void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
//     for (int vi = 0; vi < nV; vi++){
//         dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
//         Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
//         kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
//     }
// }
