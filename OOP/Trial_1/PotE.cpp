#include "PotE.hpp"
#include "PotP.hpp"
#include "CnE.hpp"
#include "CnP.hpp"
#include "mfem.hpp"
#include "MeshHandler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>
#include "mpi.h"


using namespace mfem;
using namespace std;

PotE::PotE(MeshHandler &mesh_handler, PotP &potp, CnE &cne)
    : mesh_handler(mesh_handler),
    fespace(mesh_handler.GetFESpace()),
    pmesh(mesh_handler.GetParMesh()),
    phE(fespace),
    kpl(fespace),
    RpE(fespace),
    pE0(fespace),
    cne(cne),
    potp(potp),
    B1t(fespace),
    X1v(fespace),
    B1v(fespace),
    LpCe(fespace),
    CeVn(fespace),
    CnEGridFunction(fespace),
    Dmp(cne.GetDmp())

{
    BvP = 2.9395;
    // Dmp = cne.GetDmp();
    CeVn = cne.GetCeVn();
    CnEGridFunction = cne.GetCnE();
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
    // HypreParVector LpCe(fespace);
    HypreParVector RHSl(fespace);

    Vcell = BvP - BvE;

    cout << "Vcell: " << Vcell << endl; 

}

void PotE::TimeStep(double dt){

    ParGridFunction Dmp(fespace);

    // Laplace of CnE for the RHS
    std::unique_ptr<ParBilinearForm> Kl1(new ParBilinearForm(fespace)); 
    Array<int> boundary_dofs;

    Dmp.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/DmpOOP_PotE");


    GridFunctionCoefficient cDm(&Dmp);


    Kl1->AddDomainIntegrator(new DiffusionIntegrator(cDm));
    Kl1->Assemble();
    Kl1->FormLinearSystem(boundary_dofs, phE, B1t, Kdm, X1v, B1v);

    // Vector of CnE
    CnEGridFunction.GetTrueDofs(CeVn) ;
    // CeVn.Print("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/CeVn_PotE");

    Kdm.Mult(CeVn,LpCe) ;

    // Kdm.Print("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/KdmOOP_PotE");

    // // electrolyte conductivity and RHS		
    // GridFunctionCoefficient cKe(&kpl) ;

    // std::unique_ptr<ParBilinearForm> Kl2(new ParBilinearForm(fespace)); 
            
    // Kl2->AddDomainIntegrator(new DiffusionIntegrator(cKe));
    // Kl2->Assemble();	

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

    // CGSolver cgPE(MPI_COMM_WORLD);
	// cgPE.SetRelTol(1e-7);
	// cgPE.SetMaxIter(200);		
    
    // // Solve the system using PCG with hypre's BoomerAMG preconditioner.
    // HypreBoomerAMG Mpe(Kml);
    // Mpe.SetPrintLevel(0);
    // cgPE.SetPreconditioner(Mpe);
    // cgPE.SetOperator(Kml);





}
