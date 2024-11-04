#include "Reaction.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"


Reaction::Reaction(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), AvP(*mh.GetAvP()), AvB(*mh.GetAvB())

{
    
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();


}

void Reaction::Initialize(mfem::ParGridFunction &Rx, double initial_value) {

    SetInitialReaction(Rx, initial_value);
    Rx = AvP;
    Rx *= 1.0e-8;

    // Rx.Print();


}

void Reaction::SetInitialReaction(mfem::ParGridFunction &Rx, double initial_value) {
    
    for (int i = 0; i < Rx.Size(); ++i) {
        Rx(i) = initial_value;
    }

}