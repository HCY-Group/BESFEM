#include "CnE.hpp"
#include "CnP.hpp"
#include "Reaction.hpp"
#include "Constants.hpp"
#include <iostream>
#include "MeshHandler.hpp" 
#include "mfem.hpp"
#include <iomanip>

using namespace mfem;
using namespace std;

// Reaction::Reaction(MeshHandler &mesh_handler, CnE &cne, PotE &pote)
//     : mesh_handler(mesh_handler),
//       fespace(mesh_handler.GetFESpace()),
//       Rxn(fespace),
//       Dmp(fespace),
//       kpl(fespace),
//       cne(cne),
//       CnEGridFunction(cne.GetCnE()), // defined in CnE.hpp
//       pote(pote),
//       B1t(fespace),
//       X1v(fespace),
//       B1v(fespace),
//       LpCe(fespace),
//       phE(pote.GetphE()), // defined in PotE.hpp
//       CeVn(cne.GetCeVn()), // defined in CnE.hpp
//       pmesh(mesh_handler.GetParMesh())



// {

//     pse = *mesh_handler.GetPse();
//     tc1 = pote.Gettc1();
//     tc2 = pote.Gettc2();
//     BvE = pote.GetBvE();

// }

Reaction::Reaction(MeshHandler &mesh_handler)
    : mesh_handler(mesh_handler),
      fespace(mesh_handler.GetFESpace()),
      cnp(nullptr), cne(nullptr),
      Rxn(fespace),
      CnPGridFunction(fespace),
      CnEGridFunction(fespace)



{

}

void Reaction::Initialize() {
    Rxn = 0.0;
}

void Reaction::PullCnP(CnP* cnp){
    this->cnp = cnp;
}

void Reaction::PullCnE(CnE* cne){
    this->cne = cne;
}

void Reaction::InnerLoop(){
    

		// 	// Butler-Volmer Equation for Reaction Rate
		// 	for (int vi = 0; vi < nV; vi++){
		// 		if ( AvB(vi)*dh > 0.0 ){
		// 			dPHE(vi) = phP(vi) - phE(vi);
		// 			Rxn(vi) = AvP(vi)*(Kfw(vi)*CnEGridFunction(vi)*exp(-alp*Cst1*dPHE(vi)) - \
		// 			                   Kbw(vi)*CnPGridFunction(vi)*exp( alp*Cst1*dPHE(vi)));
		// 		}
		// 	}



}

// void Reaction::TimeStep(double dt) {

    
//     nV = fespace->GetNV();

//     cout << "tc1: " << tc1 << std::endl;
//     cout << "tc2: " << tc2 << std::endl;

//     // CnEGridFunction.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/Test");

//     // electrolyte conductivity and RHS	
//     for (int vi = 0; vi < nV; vi++){
//         dffe = exp(-7.02-830*CnEGridFunction(vi)+50000*CnEGridFunction(vi)*CnEGridFunction(vi));
//         Dmp(vi) = pse(vi)*tc1*Constants::D0*dffe;
//         kpl(vi) = pse(vi)*tc2*Constants::D0*dffe*CnEGridFunction(vi);
//     }
//     GridFunctionCoefficient cDm(&Dmp);

//     // Dmp.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/DmpOOP");

//     Array<int> boundary_dofs;

//     // Laplace of CnE for the RHS
//     std::unique_ptr<ParBilinearForm> Kl1(new ParBilinearForm(fespace)); 
//     Kl1->AddDomainIntegrator(new DiffusionIntegrator(cDm));
//     Kl1->Assemble();
//     Kl1->FormLinearSystem(boundary_dofs, phE, B1t, Kdm, X1v, B1v);

//     // Vector of CnE
//     CnEGridFunction.GetTrueDofs(CeVn) ;
//     Kdm.Mult(CeVn,LpCe) ;

//     // Kdm.Print("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/KdmOOP");

//     // electrolyte conductivity and RHS		
//     GridFunctionCoefficient cKe(&kpl) ;

//     std::unique_ptr<ParBilinearForm> Kl2(new ParBilinearForm(fespace));          
//     Kl2->AddDomainIntegrator(new DiffusionIntegrator(cKe));
//     Kl2->Assemble();	

//     // assign known values to the DBC nodes	
//     ConstantCoefficient dbc_w_Coef(BvE);

//     // Dirichlet BC on the west boundary. phE
// 	Array<int> dbc_w_bdr(pmesh->bdr_attributes.Max());
// 	dbc_w_bdr = 0; dbc_w_bdr[0] = 1;
// 	// use dbc_w_bdr array to extract all node labels of Dirichlet BC
// 	Array<int> ess_tdof_list_w(0);			
// 	fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

//     phE.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 		
//     Kl2->FormLinearSystem(ess_tdof_list_w, phE, B1t, Kml, X1v, B1v);		
    
//     CGSolver cgPE(MPI_COMM_WORLD);
// 	cgPE.SetRelTol(1e-7);
// 	cgPE.SetMaxIter(200);

//     // Solve the system using PCG with hypre's BoomerAMG preconditioner.
//     HypreBoomerAMG Mpe(Kml);
//     Mpe.SetPrintLevel(0);
//     cgPE.SetPreconditioner(Mpe);
//     cgPE.SetOperator(Kml);







 
// }

// void CnE::Save() {
//     CnEGridFunction.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/OOPCnE");
// }


