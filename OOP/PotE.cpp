#include "PotE.hpp"
#include "PotP.hpp"
#include "mfem.hpp"
#include "MeshHandler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>
#include "mpi.h"


using namespace mfem;
using namespace std;

PotE::PotE(MeshHandler &mesh_handler, PotP &potp)
    : mesh_handler(mesh_handler),
    fespace(mesh_handler.GetFESpace()),
    pmesh(mesh_handler.GetParMesh()),
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


    // // assign known values to the DBC nodes	
    // ConstantCoefficient dbc_w_Coef(BvE);

    // // Dirichlet BC on the west boundary. phE
	// Array<int> dbc_w_bdr(pmesh->bdr_attributes.Max());
	// dbc_w_bdr = 0; dbc_w_bdr[0] = 1;
	// // use dbc_w_bdr array to extract all node labels of Dirichlet BC
	// Array<int> ess_tdof_list_w(0);			
	// fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);
    
    // phE.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 		
    // Kl2->FormLinearSystem(ess_tdof_list_w, phE, B1t, Kml, X1v, B1v);		
    
    // // // Solve the system using PCG with hypre's BoomerAMG preconditioner.
    // // HypreBoomerAMG Mpe(Kml);
    // // Mpe.SetPrintLevel(0);
    // // cgPE.SetPreconditioner(Mpe);
    // // cgPE.SetOperator(Kml);


}
