#include "Reaction.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"


Reaction::Reaction(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), AvP(*mh.GetAvP()), AvB(*mh.GetAvB())

{
    
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();


    Dmp = new mfem::ParGridFunction(fespace); // D_minus_plus
    kpl = new mfem::ParGridFunction(fespace); // electrolyte conductivity


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

void Reaction::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &psx1, mfem::ParGridFunction &psx2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2) {

    ElectrolyteConductivity(Cn2, psx2);
    
    mfem::GridFunctionCoefficient cDm(Dmp);
    mfem::GridFunctionCoefficient cKe(kpl);









}

void Reaction::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
    for (int vi = 0; vi < nV; vi++){

        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
        (*Dmp)(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        (*kpl)(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);

    }


}