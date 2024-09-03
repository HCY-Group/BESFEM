#include "mfem.hpp"
#include "PotP.hpp"
#include "MeshHandler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

PotP::PotP(MeshHandler &mesh_handler)
    : mesh_handler(mesh_handler),
    fespace(mesh_handler.GetFESpace()),
    phP(fespace),
    kap(fespace),
    RpP(fespace),
    pP0(fespace)
{}

void PotP::Initialize(){

    BvP = 2.9395;
    phP = BvP;

	ParBilinearForm *Kp2;

	CGSolver cgPP(MPI_COMM_WORLD);
	cgPP.SetRelTol(1e-7);
	cgPP.SetMaxIter(200);

	// force Vector
	ParLinearForm *Bp2;
	ParLinearForm Fpt(fespace);	
	HypreParVector Fpb(fespace);

	HypreParVector Xs0(fespace);
}





