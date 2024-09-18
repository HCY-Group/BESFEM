#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include "Reaction.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

Reaction::Reaction(MeshHandler &mesh_handler, Concentrations &concentrations) 
    : fespace(mesh_handler.GetFESpace()), mesh_handler(mesh_handler), 
    concentrations(concentrations)



{

    Rxn = make_unique<ParGridFunction>(fespace);


}


void Reaction::Initialize(){

    mfem::ParGridFunction &AvP = *mesh_handler.AvP;
    CreateRx(*Rxn, 0.0);
    SetAvP(*Rxn, AvP, 1.0e-08);

}

void Reaction::CreateRx(mfem::ParGridFunction &Rx, double initial_value) {

    Rx = initial_value;

}

void Reaction::SetAvP(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Av, double value) {

    Rx = Av;
    Rx *= value;
    
}
