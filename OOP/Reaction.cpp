#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include "Reaction.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

Reaction::Reaction(mfem::ParFiniteElementSpace *fe, MeshHandler &mh, Concentrations &con) 
    : fespace(fe), mesh_handler(mh), AvP(*mh.GetAvP()), concentrations(con)
    
    // , concentrations = con



{
    Rxn = new mfem::ParGridFunction(fespace);
    // concentrations = con

}


void Reaction::Initialize(){

    AvP_PGF = new ParGridFunction(fespace);	
	*AvP_PGF = AvP;

    CreateRx(*Rxn, 0.0);
    SetAvP(*Rxn, *AvP_PGF, 1.0e-08);

    // cout << "Rxn in Reactions" << endl;
    // Rxn->Print(std::cout);

}

void Reaction::CreateRx(mfem::ParGridFunction &Rx, double initial_value) {

    Rx = initial_value;

}

void Reaction::SetAvP(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Av, double value) {

    Rx = Av;
    Rx *= value;
    
}
