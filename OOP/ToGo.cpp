#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "mpi.h"

#include <cmath>
#include <memory>

using namespace std;
using namespace mfem;


int main(int argc, char *argv[])
{

   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   Hypre::Init();
   
   int myid = Mpi::WorldRank();
   
	// 1. Parse command line options.
 	const char *mesh_file = "Mesh_3x90_T3.mesh";
 	const char *dsF_file = "dsF_3x90_T3.txt"; // had "dsF_AMR_T03.txt"
	
	int order = 1;

	double dh = 0.2e-4;
	double zeta = 1.0 *0.375;				// interfacial thickness
	double eps = 1.0e-6;					// var-epsilon			
	double dt = 1.864558472553700e-01 /40;	// time step
	double tm = 0.0;						// time
	
	double t_minus = 7.619047619047619e-01; // transference number
	double D0 = 0.00489; 					// base diffusivity	
	double Frd = 96485.3365; 				// Faraday constant
	double Cst1 = 1.6021766e-19/(1.3806488e-23*300.0);
	double alp = 0.5;	
	
	double rho = 0.0501;		 	// Li site density	
	double Cr = 3.0;				// C-rate
	double Vsr = 0.2;				// voltage scanning rate
	double Vcut = 2.7; 				// cut-off voltage
	


	// // =====================================
	// // 	  __  __           _     
	// // 	 |  \/  |         | |    
	// // 	 | \  / | ___  ___| |__  
	// // 	 | |\/| |/ _ \/ __| '_ \ 		//
	// // 	 | |  | |  __/\__ \ | | |
	// // 	 |_|  |_|\___||___/_| |_|
	// // =====================================                        
                         
	
	// Read the mesh from the given mesh file.
	Mesh gmesh(mesh_file);
	gmesh.EnsureNCMesh(true);	

	// Create global FE space for distance function.
	H1_FECollection gFec(order, gmesh.Dimension());	
	FiniteElementSpace gFespace(&gmesh, &gFec);
	
	// Read global distance function
	GridFunction gDsF(&gFespace);   
	int Onm = gDsF.Size(); 
	ifstream myfile;
	myfile.open(dsF_file);
	for(int gi = 0; gi < Onm ; gi++){
		myfile >> gDsF(gi); 
	}   
	myfile.close(); 	

	// west boundary size
	Vector Rmin, Rmax;
	gmesh.GetBoundingBox(Rmin,Rmax);
	double L_w = Rmax(1) - Rmin(1);
	
	// local (parallel) mesh 
	ParMesh pmesh(MPI_COMM_WORLD, gmesh);
	pmesh.Save("/mnt/scratch/brandlan/3x90/Pmesh_AMR_T4");
	
	int nV = pmesh.GetNV();					// number of vertices
	int nE = pmesh.GetNE();					// number of elements
	int nC = pow(2,pmesh.Dimension());		// number of corner vertices
	
	Array<double> VtxVal(nC) ;				// values of corner vertices
	
	// element area
	Vector EVol(nE);
	for (int ei = 0; ei < nE; ei++){
		EVol(ei) = pmesh.GetElementVolume(ei);	
	} 	

	cout << "Number of vertices: " << nV << endl;
    cout << "Number of elements: " << nE << endl;
	cout << "West boundary size: " << L_w << endl;
	
	
	// 	=====================================================================
	// 	//   _____                        _                               
	// 	//  |  __ \                      (_)                              
	// 	//  | |  | | ___  _ __ ___   __ _ _ _ __    _ __   __ _ _ __ __ _ 
	// 	//  | |  | |/ _ \| '_ ` _ \ / _` | | '_ \  | '_ \ / _` | '__/ _` |
	// 	//  | |__| | (_) | | | | | | (_| | | | | | | |_) | (_| | | | (_| |
	// 	//  |_____/ \___/|_| |_| |_|\__,_|_|_| |_| | .__/ \__,_|_|  \__,_|
	// 	//                                         | |                    
	// 	//                                         |_|                    
	// 	=====================================================================
	

	// Create local FE space.
	H1_FECollection fec(order, pmesh.Dimension());	
	ParFiniteElementSpace fespace(&pmesh, &fec);
	// 	HYPRE_BigInt total_num_dofs = fespace.GlobalTrueVSize();		
	
	// Map local to global element indices.
	Array<HYPRE_BigInt> E_L2G;    
	pmesh.GetGlobalElementIndices(E_L2G);		
		
	Array<int> gVTX(nC);					// global indices of corner vertices				
	Array<int> VTX(nC);						// local indices of corner vertices
	int gei;								// global element indices

	// Local (parallel) GridFunction
	ParGridFunction dsF(&fespace);	

	// Map local distance function from global one
	for (int ei = 0; ei < nE; ei++){
		gei = E_L2G[ei];
	// 	gei = pmesh.GetGlobalElementNum(ei);
		
		gmesh.GetElementVertices(gei,gVTX);
		pmesh.GetElementVertices(ei,VTX);
	
		for (int vi = 0; vi < nC; vi++){
			dsF(VTX[vi]) = gDsF(gVTX[vi]);
		}			
	}	
	
	// Local domain parameters
	ParGridFunction psi(&fespace);
	ParGridFunction pse(&fespace);	
	ParGridFunction AvP(&fespace);	

	// interpolate domain parameter from distance function
    for(int vi = 0; vi < nV; vi++){
        psi(vi) = 0.5*(1.0 + tanh(dsF(vi)/(zeta*dh)));  
        pse(vi) = 1.0 - psi(vi);
        AvP(vi) = -(pow(tanh(dsF(vi)/(zeta*dh)),2) - 1.0)/(2*zeta*dh);
        
        if (psi(vi) < eps){psi(vi) = eps;}   
        if (pse(vi) < eps){pse(vi) = eps;}         
    } 	
	
	ParGridFunction AvB(&fespace);
	AvB = AvP;
	for (int vi = 0; vi < nV; vi++){
		 // remove small values of AvP  
		if (AvP(vi)*dh < 1.0e-3){AvP(vi) = 0.0;} 
		if (AvB(vi)*dh < 1.0e-6){AvB(vi) = 0.0;}  		 
	}	
	
	// Total psi throughout the global domain
	double tPsi = 0.0;
	double val = 0.0;
	
	
	Vector EAvg(nE);						// element average	
	for (int ei = 0; ei < nE; ei++){
		psi.GetNodalValues(ei,VtxVal) ;
		val = 0.0;
		for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
		EAvg(ei) = val/nC;
		tPsi += EAvg(ei)*EVol(ei) ;
	} 	
	
	
	double gtPsi;
	MPI_Allreduce(&tPsi, &gtPsi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
	
	// Total pse throughout the global domain	
	double tPse = 0.0;	
	for (int ei = 0; ei < nE; ei++){
		pse.GetNodalValues(ei,VtxVal) ;
		val = 0.0;
		for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
		EAvg(ei) = val/nC;
		tPse += EAvg(ei)*EVol(ei) ;
	} 	
	
	double gtPse;
	MPI_Allreduce(&tPse, &gtPse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Target Current
	double trgI = tPsi*rho*(0.9-0.3)/(3600.0/Cr);
	double gTrgI = 0.0;
	MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	double sCrnt = 0.0;	
	double gCrnt = 0.0;	
		
	// Some containers that will be used later.
	Array<int> boundary_dofs;					// nature boundary
	
	// Neumann BC on the west boundary. CnE
	Array<int> nbc_w_bdr(pmesh.bdr_attributes.Max());
	nbc_w_bdr = 0; nbc_w_bdr[0] = 1;

	// Dirichlet BC on the east boundary. phP
	Array<int> dbc_e_bdr(pmesh.bdr_attributes.Max());
	dbc_e_bdr = 0; dbc_e_bdr[2] = 1;
	// use dbc_e_bdr array to extract all node labels of Dirichlet BC
	Array<int> ess_tdof_list_e(0);			
	fespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);
	
	// Dirichlet BC on the west boundary. phE
	Array<int> dbc_w_bdr(pmesh.bdr_attributes.Max());
	dbc_w_bdr = 0; dbc_w_bdr[0] = 1;
	// use dbc_w_bdr array to extract all node labels of Dirichlet BC
	Array<int> ess_tdof_list_w(0);			
	fespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);	
		
	// delete global mesh because it is longer used
	gmesh.Clear();
		
	MPI_Barrier(MPI_COMM_WORLD);

	cout << "Total Psi: " << gtPsi << endl;
    cout << "Total Pse: " << gtPse << endl;
    cout << "Target Current: " << gTrgI << endl;
	
		
// 	// // 	=============================
// 	// // 	   _____       _____  
// 	// // 	  / ____|     |  __ \ 		//
// 	// // 	 | |     _ __ | |__) |
// 	// // 	 | |    | '_ \|  ___/ 
// 	// // 	 | |____| | | | |     
// 	// // 	  \_____|_| |_|_|     
// 	// // 	==============================

	// initial condition
	ParGridFunction CnP(&fespace);
	double Cp0 = 0.3;					// initial value
	CnP = Cp0;	
	
	// degree of lithiation
	double Xfr = 0.0;

	ParGridFunction TmpF(&fespace);
	TmpF = CnP;
	TmpF *= psi;
	double lSum = 0.0;
	for (int ei = 0; ei < nE; ei++){	  
		TmpF.GetNodalValues(ei,VtxVal);
		val = 0.0;
		for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
		EAvg(ei) = val/nC;		 		  
		lSum += EAvg(ei)*EVol(ei);	
	} 	


    cout << "Sum: " << lSum << std::endl;


	double gSum;
	MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
	Xfr = gSum/gtPsi;	
		 

	// SBM mass matrix
	HypreParMatrix Mmatp;
	std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(&fespace));
 	GridFunctionCoefficient cPs(&psi) ;
 	Mt->AddDomainIntegrator(new MassIntegrator(cPs)); 	
 	Mt->Assemble();
 	Mt->FormSystemMatrix(boundary_dofs, Mmatp);
	
	HypreSmoother Mp_prec;
	CGSolver Mp_solver(MPI_COMM_WORLD);	

	Mp_solver.iterative_mode = false;
	Mp_solver.SetRelTol(1e-7);
	Mp_solver.SetAbsTol(0);
	Mp_solver.SetMaxIter(200);
	Mp_solver.SetPrintLevel(0);
	Mp_prec.SetType(HypreSmoother::Jacobi);
	Mp_solver.SetPreconditioner(Mp_prec);
	Mp_solver.SetOperator(Mmatp);

	HypreParMatrix *Tmatp;

	// SBM stiffness matrix
	ParGridFunction Dp(&fespace);		

	HypreParMatrix Kmatp;
	//ParBilinearForm *Kc2;	
	
	// force vector
	ParGridFunction Rxc(&fespace);
	ParLinearForm *Bc2;
	ParLinearForm Fct(&fespace);
	HypreParVector Fcb(&fespace);	
	
 	// create a Vector for CnP
	HypreParVector CpV0(&fespace), CpVn(&fespace), RHCp(&fespace);		
	
	int nDof = CpV0.Size();	
	
	// Vector of psi
	HypreParVector PsVc(&fespace);
	psi.GetTrueDofs(PsVc);	

	cout << "Degree of Lithiation: " << Xfr << std::endl;	


	ParGridFunction Rxn(&fespace);
	Rxn = 0.0;
	Rxn = AvP;
	Rxn *= 1.0e-10;	
	
	//containers used later.
	HypreParVector X1v(&fespace);

	string stri;

		  
		
int t = 0;
for (int t = 0; t < 10 + 1; t++){
	
	
		// // 	=============================
		// // 	   _____       _____  
		// // 	  / ____|     |  __ \ 		//
		// // 	 | |     _ __ | |__) |
		// // 	 | |    | '_ \|  ___/ 
		// // 	 | |____| | | | |     
		// // 	  \_____|_| |_|_|     
		// // 	==============================
	
			
		Rxc = Rxn;
		Rxc /= rho;	
		GridFunctionCoefficient cAp(&Rxc) ;	

		// force term	
		std::unique_ptr<ParLinearForm> Bc2(new ParLinearForm(&fespace));
		Bc2->AddDomainIntegrator(new DomainLFIntegrator(cAp));
		Bc2->Assemble();
		// Move the contents of Bc2 into Fct
		Fct = std::move(*Bc2);
		
		// diffusivity in the particles
		for (int vi = 0; vi < nV; vi++ ){
			// appendix equation A-19
			Dp(vi) = psi(vi)*(0.0277-0.084*CnP(vi) + 0.1003*CnP(vi)*CnP(vi))*1.0e-8;
			if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}
		}	
		GridFunctionCoefficient cDp(&Dp) ;			
		
		// K matrix		
		std::unique_ptr<ParBilinearForm> Kc2(new ParBilinearForm(&fespace)); 
		HypreParMatrix Kmatp; 

   		Kc2->AddDomainIntegrator(new DiffusionIntegrator(cDp));
   		Kc2->Assemble();
   		Kc2->FormLinearSystem(boundary_dofs, CnP, Fct, Kmatp, X1v, Fcb);
   		Fcb *= dt;

		// T matrix
		Tmatp = Add(1.0, Mmatp, -dt, Kmatp);	
		
		// vector of CnP				
		CnP.GetTrueDofs(CpV0);
		
		Tmatp->Mult(CpV0,RHCp);
		RHCp += Fcb;		

		// time stepping
		Mp_solver.Mult(RHCp,CpVn) ;		

		// Update only the solid region
		for (int p = 0; p < nDof; p++){
			if (PsVc(p) < 1.0e-5){
				CpVn(p) = Cp0;}
		}
		
		// Recover the GridFunction from Vector.
		CnP.Distribute(CpVn);

		// Degree of lithiation
		TmpF = CnP;
		TmpF *= psi;
		lSum = 0.0;
		for (int ei = 0; ei < nE; ei++){	  
			TmpF.GetNodalValues(ei,VtxVal);
			val = 0.0;
			for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
			EAvg(ei) = val/nC;		 		  
			lSum += EAvg(ei)*EVol(ei);	
		} 	
		MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
		Xfr = gSum/gtPsi;	

		delete Tmatp;
		
 	}

	
	CnP.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/CnP");
		
	
   return 0;
}	