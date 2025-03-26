/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "PotE.hpp"
#include "../code/Constants.hpp"
#include "mfem.hpp"
#include "CnE.hpp"

double BvE = 0.0; ///< Global variable for the boundary value of electrolyte potential

PotE::PotE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_w_bdr(geo.dbc_w_bdr), gtPse(para.gtPse), ess_tdof_list_w(geo.ess_tdof_list_w)
    
    {
    cgPE_solver = std::make_unique<mfem::CGSolver>(MPI_COMM_WORLD);
    RpE = std::make_unique<mfem::ParGridFunction>(fespace.get()); // reaction rate for particle;
    Flb = mfem::HypreParVector(fespace.get());
    Dmp = std::make_unique<mfem::ParGridFunction>(fespace.get()); // D_minus_plus
    kpl = std::make_unique<mfem::ParGridFunction>(fespace.get()); // electrolyte conductivity
    B1t = mfem::ParLinearForm(fespace.get());
    X1v = mfem::HypreParVector(fespace.get());
    B1v = mfem::HypreParVector(fespace.get());
    RHSl = mfem::HypreParVector(fespace.get());
    LpCe = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    CeVn = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    }

void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value)
{
    BvE = initial_value; // Set the boundary value.
    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials
    Potentials::SetUpSolver(*cgPE_solver, 1e-7, 80); // Configure the solver
    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);
}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{
    ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity
    auto cDm = std::make_shared<mfem::GridFunctionCoefficient>(Dmp.get());
    auto cKe = std::make_shared<mfem::GridFunctionCoefficient>(kpl.get());

    Solver::StiffnessMatrix(cDm, boundary_dofs, potential, B1t, Kdm, B1v);
    Cn.GetTrueDofs(*CeVn); 
    Kdm->Mult(*CeVn, *LpCe); // Multiply concentration by diffusivity matrix
    
    Potentials::ImplementBoundaryConditions(dbc_w_Coef, BvE, potential, dbc_w_bdr); // Apply boundary conditions
    Solver::StiffnessMatrix(cKe, ess_tdof_list_w, potential, B1t, KmE, B1v);
    Potentials::PCG_Solver(Mpe, *cgPE_solver, *KmE); // Solve the system
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

    // Create dummy values to use Force Term Function 
    static mfem::Array<int> dummy_boundary;
    static mfem::ConstantCoefficient coef1(0.0);
    static mfem::ConstantCoefficient coef2(0.0);
    static mfem::ProductCoefficient dummy_coef(coef1, coef2);

    // Assemble the force term without applying boundary conditions
    Solver::ForceTerm(*RpE, ftPotE, dummy_boundary, dummy_coef, false); // false since not applying BCs
    
    Potentials::ForceVector(*K, ess_tdof_list_w, phx, ftPotE, *KmE, X1v, Flb, dbc_w_Coef, dbc_w_bdr); // Force vector
    RHSl = Flb;
    RHSl += *LpCe;
    Potentials::ErrorCalculation(phx, *cgPE_solver, RHSl, psx, error_E, gerror, gtPse); // Compute global error
}