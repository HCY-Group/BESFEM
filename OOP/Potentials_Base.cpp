#include "Potentials_Base.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"


Potentials::Potentials(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh)

{
    
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();




}


void Potentials::SetInitialPotentials(mfem::ParGridFunction &ph, double initial_value) {
    
    for (int i = 0; i < ph.Size(); ++i) {
        ph(i) = initial_value;
    }

}

void Potentials::SetUpSolver(mfem::CGSolver &solver, double value_1, double value_2) {
    
    solver.SetRelTol(value_1);
    solver.SetMaxIter(value_2);

}


void Potentials::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {

    Rx2 = Rx1;
    Rx2 *= value;

}

void Potentials::ForceTerm(mfem::ParGridFunction &Rx2, mfem::ParLinearForm &Fxx) {

    std::unique_ptr<ParLinearForm> Bx2(new ParLinearForm(fespace));	

    Rxx = new ParGridFunction(fespace);
    *Rxx = Rx2;
    cXx = new GridFunctionCoefficient(Rxx);

    Bx2->AddDomainIntegrator(new DomainLFIntegrator(*cXx));
    Bx2->Assemble();
    Fxx = std::move(*Bx2);

    Bx2->Print(std::cout);

}