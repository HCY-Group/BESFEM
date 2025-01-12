/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "PotE.hpp"
#include "mfem.hpp"
#include "CnE.hpp"

double BvE = 0.0; ///< Global variable for the boundary value of electrolyte potential

PotE::PotE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh), fespace(fe), dbc_w_bdr(mh.dbc_w_bdr), gtPse(mh.gtPse), ess_tdof_list_w(mh.ess_tdof_list_w)
    
    {

    cgPE_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    RpE = new ParGridFunction(fespace); // reaction rate for particle;

    Flb = mfem::HypreParVector(fespace);

    Dmp = new mfem::ParGridFunction(fespace); // D_minus_plus
    kpl = new mfem::ParGridFunction(fespace); // electrolyte conductivity

    Kl1 = std::make_unique<mfem::ParBilinearForm>(fespace); // Use make_unique
    Kl2 = std::make_unique<mfem::ParBilinearForm>(fespace); // Use make_unique

    B1t = mfem::ParLinearForm(fespace);
    X1v = mfem::HypreParVector(fespace);
    B1v = mfem::HypreParVector(fespace);
    RHSl = mfem::HypreParVector(fespace);

    LpCe = new mfem::HypreParVector(fespace);
    CeVn = new mfem::HypreParVector(fespace);

    // std::cout << "fespace address in PotE: " << fespace << std::endl;


    }

mfem::CGSolver* PotE::cgPE_solver = nullptr; // static variable to be used in reaction


void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    BvE = initial_value; // Set the boundary value.

    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials
    Potentials::SetUpSolver(*cgPE_solver, 1e-7, 80); // Configure the solver

    // std::cout << "fespace address in PotE_init: " << fespace << std::endl;


    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

    
}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx){

    ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    mfem::GridFunctionCoefficient cDm(Dmp); 
    mfem::GridFunctionCoefficient cKe(kpl); 
    
    Kl1 = std::make_unique<mfem::ParBilinearForm>(fespace);  // Initialize as a member variable

    // Kl1->Update();

    Potentials::KMatrix(*Kl1, cDm, boundary_dofs, phx, B1t, Kdm, X1v, B1v); // Diffusivity matrix

    Cn.GetTrueDofs(*CeVn); 
    Kdm.Mult(*CeVn, *LpCe); // Multiply concentration by diffusivity matrix

    // Kl2->Update();

    Kl2 = std::make_unique<mfem::ParBilinearForm>(fespace);  // Initialize as a member variable

    Potentials::ImplementBoundaryConditions(dbc_w_Coef, BvE, phx, dbc_w_bdr); // Apply boundary conditions
    Potentials::KMatrix(*Kl2, cKe, ess_tdof_list_w, phx, B1t, KmE, X1v, B1v); // Conductivity matrix
    Potentials::PCG_Solver(Mpe, *cgPE_solver, KmE); // Solve the system

}

void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
    for (int vi = 0; vi < nV; vi++){

        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        (*Dmp)(vi) = psx(vi) * tc1 * Constants::D0 * dffe; // Update diffusivity
        (*kpl)(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi); // Update conductivity

    }

}

void PotE::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) 

{

    Potentials::CreateReaction(Rx, *RpE, -1.0); // Create reaction field

    Potentials::ForceTerm(*RpE, ftPotE);  // Assemble the force term
    Potentials::ForceVector(*Kl2, ess_tdof_list_w, phx, ftPotE, KmE, X1v, Flb, dbc_w_Coef, dbc_w_bdr); // Force vector

    RHSl = Flb;
    RHSl += *LpCe;


    Potentials::ErrorCalculation(phx, *cgPE_solver, RHSl, psx, error_E, gerror, gtPse); // Compute global error

}