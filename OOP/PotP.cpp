/**
 * @file PotP.cpp
 * @brief Implementation of the potential class for particle potential simulations.
 */

#include "PotP.hpp"
#include "mfem.hpp"

double BvP = 0.0; ///< Global variable for the boundary value of particle potential

PotP::PotP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh), fespace(fe), dbc_e_bdr(mh.dbc_e_bdr), gtPsi(mh.gtPsi), ess_tdof_list_e(mh.ess_tdof_list_e)
    
    {

    cgPP_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    RpP = new ParGridFunction(fespace); // Initialize reaction rate grid function

    Fpb = mfem::HypreParVector(fespace); // Initialize force term vector
    kap = new mfem::ParGridFunction(fespace); // Initialize conductivity field
    Kp2 = std::make_unique<mfem::ParBilinearForm>(fespace); // Initialize bilinear form

    B1t = mfem::ParLinearForm(fespace);
    X1v = mfem::HypreParVector(fespace);
    B1v = mfem::HypreParVector(fespace);

    }

mfem::CGSolver* PotP::cgPP_solver = nullptr; // static variable to be used in reaction


void PotP::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    
    BvP = initial_value; // Set the boundary value
    
    Potentials::SetInitialPotentials(ph, BvP); // Initialize potentials
    Potentials::SetUpSolver(*cgPP_solver, 1e-7, 82); // Configure the solver

    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e); // Identify essential DoFs


    
}


void PotP::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx)


{
    ParticleConductivity(Cn, psx); // Update conductivity

    mfem::GridFunctionCoefficient cKp(kap); // Wrap conductivity as a coefficient

    Potentials::ImplementBoundaryConditions(dbc_e_Coef, BvP, phx, dbc_e_bdr); // Apply BCs

    Kp2 = std::make_unique<ParBilinearForm>(fespace);  // Initialize as a member variable

    // Kp2->Update();

    Potentials::KMatrix(*Kp2, cKp, ess_tdof_list_e, phx, B1t, KmP, X1v, B1v); // Assemble system matrix
    Potentials::PCG_Solver(Mpp, *cgPP_solver, KmP); // Solve the system

}



void PotP::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){

    for (int vi = 0; vi < nV; vi++){

        (*kap)(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);

    }

}

void PotP::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) 

{

    Potentials::CreateReaction(Rx, *RpP, Constants::Frd); // Create reaction field
    Potentials::ForceTerm(*RpP, ftPotP); // Assemble the force term
    Potentials::ForceVector(*Kp2, ess_tdof_list_e, phx, ftPotP, KmP, X1v, Fpb, dbc_e_Coef, dbc_e_bdr); // Force vector

    Potentials::ErrorCalculation(phx, *cgPP_solver, Fpb, psx, error_P, gerror, gtPsi); // Compute global error


}


