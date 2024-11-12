#include "PotP.hpp"
#include "mfem.hpp"


PotP::PotP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh)
    
    {

    cgPP_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    RpP = new ParGridFunction(fespace); // reaction rate for particle;

    Fpb = mfem::HypreParVector(fespace);

    }

mfem::CGSolver* PotP::cgPP_solver = nullptr; // static variable to be used in reaction


void PotP::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    Potentials::SetInitialPotentials(ph, initial_value);
    Potentials::SetUpSolver(*cgPP_solver, 1e-7, 82);

    
}

void PotP::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx) 

{

    Potentials::CreateReaction(Rx, *RpP, Constants::Frd);
    Potentials::ForceTerm(*RpP, ftPotP);





}