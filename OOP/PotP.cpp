#include "PotP.hpp"
#include "mfem.hpp"


PotP::PotP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh)
    
    {

    cgPP_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    RpP = new ParGridFunction(fespace); // reaction rate for particle;

    Fpb = mfem::HypreParVector(fespace);
    kap = new mfem::ParGridFunction(fespace); // particle conductivity
    Kp2 = new mfem::ParBilinearForm(fespace);

    B1t = new mfem::ParLinearForm(fespace);
    X1v = new mfem::HypreParVector(fespace);
    B1v = new mfem::HypreParVector(fespace);


    }

mfem::CGSolver* PotP::cgPP_solver = nullptr; // static variable to be used in reaction


void PotP::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    Potentials::SetInitialPotentials(ph, initial_value);
    Potentials::SetUpSolver(*cgPP_solver, 1e-7, 82);

    
}


void PotP::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx)


{
    ParticleConductivity(Cn, psx);

    mfem::GridFunctionCoefficient cKp(kap);

    mfem::Array<int> ess_tdof_list_e = mesh_handler.ess_tdof_list_e;

    //Potentials::ImplementBoundaryConditions();
    Potentials::KMatrix(*Kp2, cKp, ess_tdof_list_e, phx, *B1t, KmP, *X1v, *B1v);
    Potentials::PCG_Solver(Mpp, *cgPP_solver, KmP);

}



void PotP::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){

    for (int vi = 0; vi < nV; vi++){

        (*kap)(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);

    }

}

void PotP::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx) 

{

    Potentials::CreateReaction(Rx, *RpP, Constants::Frd);
    Potentials::ForceTerm(*RpP, ftPotP);
    // Potentials::ForceVector(Kp2, ess_tdof_list_e, phx, ftPotP, KmP, X1v, Fpb, dbc_e_Coef, dbc_e_bdr);





}