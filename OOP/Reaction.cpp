#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include "Reaction.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

Reaction::Reaction(mfem::ParFiniteElementSpace *fe, MeshHandler &mesh_handler, Concentrations &concentrations) 
    : fespace(fe), mesh_handler(mesh_handler), concentrations(concentrations), AvP(*mesh_handler.GetAvP())



{
    Rxn = new ParGridFunction(fespace);

}


void Reaction::Initialize(){

    // mfem::ParGridFunction &AvP = *mesh_handler.AvP;

    AvP_PGF = new ParGridFunction(fespace);	
	*AvP_PGF = AvP;

    CreateRx(*Rxn, 0.0);
    SetAvP(*Rxn, *AvP_PGF, 1.0e-08);

    // Rxn->Print(std::cout);

}

void Reaction::CreateRx(mfem::ParGridFunction &Rx, double initial_value) {

    Rx = initial_value;

}

void Reaction::SetAvP(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Av, double value) {

    Rx = Av;
    Rx *= value;
    
}
