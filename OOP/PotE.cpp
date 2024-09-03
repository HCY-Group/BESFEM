#include "PotE.hpp"
#include "PotP.hpp"
#include "mfem.hpp"
#include "MeshHandler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>


using namespace mfem;
using namespace std;

PotE::PotE(MeshHandler &mesh_handler, PotP &potp)
    : mesh_handler(mesh_handler),
    fespace(mesh_handler.GetFESpace()),
    phE(fespace),
    kpl(fespace),
    RpE(fespace),
    pE0(fespace),
    // Dmp(fespace), // not used in PotE
    potp(potp)

{
    BvP = potp.GetBvP();
}

void PotE::Initialize(){

    tc1 = (2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
    tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    BvE = -1.0;

    phE = BvE;

	ParBilinearForm *Kl1;
    ParBilinearForm *Kl2;

	CGSolver cgPE(MPI_COMM_WORLD);
	cgPE.SetRelTol(1e-7);
	cgPE.SetMaxIter(200);

	// force Vector
	ParLinearForm *Bl2;
	ParLinearForm Flt(fespace);	
	HypreParVector Flb(fespace);

	HypreParVector Xe0(fespace);
    HypreParVector LpCe(fespace);
    HypreParVector RHSl(fespace);

    Vcell = BvP - BvE;

    // cout << "BvP in PotE: " << BvP << endl; 
    // cout << "BvE: " << BvE << endl; 

    cout << "Vcell: " << Vcell << endl; 
}















// void InitializePotE(ParFiniteElementSpace &fespace)
// {
//     double t_minus = 7.619047619047619e-01; // transference number
//     double Cst1 = 1.6021766e-19/(1.3806488e-23*300.0);
//     double tc1 = (2*t_minus - 1.0) / (2 * t_minus * (1.0 - t_minus));
//     double tc2 = 1.0 / (2 * t_minus * (1.0 - t_minus)) * Cst1;
//     double dffe;

//     ParGridFunction phE(&fespace);    // electropotential in electrolyte
//     ParGridFunction Dmp(&fespace);    // D_minus_plus
//     ParGridFunction kpl(&fespace);    // electrolyte conductivity
//     ParGridFunction RpE(&fespace);    // reaction rate for electrolyte
//     ParGridFunction pE0(&fespace);

//     BvE = -1.0;
//     phE = BvE;

//     // stiffness matrix
//     HypreParMatrix Kml;
//     ParBilinearForm *Kl2;

//     // force vector
//     ParLinearForm *Bl2;
//     ParLinearForm Flt(&fespace);
//     HypreParVector Flb(&fespace);

//     // Laplace matrix
//     HypreParMatrix Kdm;
//     ParBilinearForm *Kl1;

//     CGSolver cgPE(MPI_COMM_WORLD);
//     cgPE.SetRelTol(1e-7);
//     cgPE.SetMaxIter(200);

//     //cout << "BvETEST: " << BvE << std::endl;

// }
