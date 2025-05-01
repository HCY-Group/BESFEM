/**
 * @file PotP.cpp
 * @brief Implementation of the potential class for particle potential simulations.
 */

#include "PotP.hpp"
#include "../code/Constants.hpp"
#include "mfem.hpp"
#include <optional>


// double BvP = 0.0; ///< Global variable for the boundary value of particle potential
double BvP = -0.1; ///< Global variable for the boundary value of particle potential

PotP::PotP(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_e_bdr(geo.dbc_e_bdr), gtPsi(para.gtPsi), ess_tdof_list_e(geo.ess_tdof_list_e)
    
    {
    cgPP_solver = std::make_unique<mfem::CGSolver>(MPI_COMM_WORLD);
    RpP = std::make_unique<mfem::ParGridFunction>(fespace.get()); // Initialize reaction rate grid function
    Fpb = mfem::HypreParVector(fespace.get()); // Initialize force term vector
    kap = std::make_unique<mfem::ParGridFunction>(fespace.get()); // Initialize conductivity field
    B1t = mfem::ParLinearForm(fespace.get());
    X1v = mfem::HypreParVector(fespace.get());
    B1v = mfem::HypreParVector(fespace.get());
    }


void PotP::Initialize(mfem::ParGridFunction &ph, double initial_value)
{
    BvP = initial_value; // Set the boundary value
    Potentials::SetInitialPotentials(ph, BvP); // Initialize potentials
    Potentials::SetUpSolver(*cgPP_solver, 1e-7, 82); // Configure the solver
    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e); // Identify essential DoFs
}


void PotP::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{
    ParticleConductivity(Cn, psx); // Update conductivity
    auto cKp = std::make_shared<mfem::GridFunctionCoefficient>(kap.get());

    Potentials::ImplementBoundaryConditions(dbc_e_Coef, BvP, potential, dbc_e_bdr); // Apply BCs
    SolverSteps::StiffnessMatrix(cKp, ess_tdof_list_e, potential, B1t, KmP, B1v);
    Potentials::PCG_Solver(Mpp, *cgPP_solver, *KmP); // Solve the system
}

void PotP::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){
    for (int vi = 0; vi < nV; vi++){
        // (*kap)(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);
        (*kap)(vi) = 3.3 * psx(vi);
    }
}

void PotP::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) 
{
    Potentials::CreateReaction(Rx, *RpP, Constants::Frd); // Create reaction field

    // Assemble the force term without applying boundary conditions
    SolverSteps::ForceTerm(*RpP, ftPotP); // false since not applying BCs

    Potentials::ForceVector(*K, ess_tdof_list_e, phx, ftPotP, *KmP, X1v, Fpb, dbc_e_Coef, dbc_e_bdr); // Force vector
    Potentials::ErrorCalculation(phx, *cgPP_solver, Fpb, psx, error_P, gerror, gtPsi); // Compute global error
}


