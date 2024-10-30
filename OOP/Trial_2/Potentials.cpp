#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include "Potentials.hpp"
#include <fstream>
#include <iostream>

#include <HYPRE.h>
#include <HYPRE_utilities.h>
#include "mfem.hpp"


Potentials::Potentials(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
: mesh_handler(mh), fespace(fe), pmesh(pm)

{

    // nE = mesh_handler.GetNE();
    // nC = mesh_handler.GetNC();
    // nV = mesh_handler.GetNV();
    
    phP = new mfem::ParGridFunction(fespace); // electropotential in particle
    phE = new mfem::ParGridFunction(fespace); // electropotential in electrolyte

    cgPP_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    cgPE_solver = new mfem::CGSolver(MPI_COMM_WORLD);


    // kap = new ParGridFunction(fespace); // conductivity in particle
    // RpP = new ParGridFunction(fespace); // reaction
    // pP0 = new ParGridFunction(fespace); // values before iteration

}

void Potentials::InitializePotP() {

    BvP = 2.9395;
    CreatePotentials(*phP, BvP);
    Solver(*cgPP_solver, 1e-7, 82); 

}

void Potentials::InitializePotE() {

    BvE = -1.0;
    CreatePotentials(*phE, BvE);
    Solver(*cgPE_solver, 1e-7, 80);

    Vcell = BvP - BvE;

}


void Potentials::CreatePotentials(mfem::ParGridFunction &ph, double initial_value){

    for (int i = 0; i < ph.Size(); ++i) {
        ph(i) = initial_value;  // Set all values of Cn to initial_value
    }  

}

void Potentials::Solver(mfem::CGSolver &solver, double value_1, double value_2) {

    solver.SetRelTol(value_1);
    solver.SetMaxIter(value_2);

}

void Potentials::SetReaction(Reaction *reaction_) {
    reaction = reaction_;  // Store the Reaction object for later use
}

