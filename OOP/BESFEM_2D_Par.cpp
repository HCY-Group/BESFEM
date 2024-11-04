#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "mpi.h"

#include <cmath>

// #include <stdio.h>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{

   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   Hypre::Init();
   
   int myid = Mpi::WorldRank();
   
	// 1. Parse command line options.
//	string mesh_file = "Mesh/Smesh_80x90_AMR_D01.mesh";
//	string dsF_file =  "Mesh/dsF_81x91_AMR_D01.txt";
		
 	const char *mesh_file = "Mesh_3x90_T3.mesh";
 	const char *dsF_file = "dsF_3x90_T3.txt"; // had "dsF_AMR_T03.txt" 


// 	bool visualization = true;
	
	int order = 1;

	double dh = 0.5e-4;
	double zeta = 1.0 *0.375;				// interfacial thickness
	double thres = 1.0e-3;					// AvP musk
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
	pmesh.Save("Output/Pmesh_AMR_D01");
	
	int nV = pmesh.GetNV();					// number of vertices
	int nE = pmesh.GetNE();					// number of elements
	int nC = pow(2,pmesh.Dimension());		// number of corner vertices
	
	Array<double> VtxVal(nC) ;				// values of corner vertices
	
	// element area
	Vector EVol(nE);
	for (int ei = 0; ei < nE; ei++){
		EVol(ei) = pmesh.GetElementVolume(ei);	
	} 	
	
	
// 	// 	=====================================================================
// 	// 	//   _____                        _                               
// 	// 	//  |  __ \                      (_)                              
// 	// 	//  | |  | | ___  _ __ ___   __ _ _ _ __    _ __   __ _ _ __ __ _ 
// 	// 	//  | |  | |/ _ \| '_ ` _ \ / _` | | '_ \  | '_ \ / _` | '__/ _` |
// 	// 	//  | |__| | (_) | | | | | | (_| | | | | | | |_) | (_| | | | (_| |
// 	// 	//  |_____/ \___/|_| |_| |_|\__,_|_|_| |_| | .__/ \__,_|_|  \__,_|
// 	// 	//                                         | |                    
// 	// 	//                                         |_|                    
// 	// 	=====================================================================
	

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
// 		gei = pmesh.GetGlobalElementNum(ei);
		
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
	
	// for (int i = 0; i < psi.Size(); ++i) {
    // cout << "psi[" << i << "] = " << psi(i) << std::endl;
    // }	
	
	ParGridFunction AvB(&fespace);
	AvB = AvP;
	for (int vi = 0; vi < nV; vi++){
		 // remove small values of AvP  
		if (AvP(vi)*dh < thres ){AvP(vi) = 0.0;} 
		if (AvB(vi)*dh < 1.0e-6){AvB(vi) = 0.0;}  		 
	}
 	
	// std::cout << "AvB: " << AvB << std::endl;
	
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
	
// 	double sCrnt = 0.0;	
// 	double gCrnt = 0.0;	
		
	// Some containers that will be used later.
	Array<int> boundary_dofs;					// nature boundary
	
	// // Neumann BC on the west boundary. CnE
	Array<int> nbc_w_bdr(pmesh.bdr_attributes.Max());
	nbc_w_bdr = 0; 
	nbc_w_bdr[0] = 1;

    // // Printing the values
    // std::cout << "nbc_w_bdr values: ";
    // for (int i = 0; i < nbc_w_bdr.Size(); i++) {
    //     std::cout << nbc_w_bdr[i] << " ";
    // }
    // std::cout << std::endl;

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
	
// // 	ConstantCoefficient one(1.0);		
		
	
// 	// delete global mesh because it is longer used
// 	gmesh.Clear();
		
	MPI_Barrier(MPI_COMM_WORLD);
	
		
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

	// CnP.Print(std::cout);
	
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

	double gSum;
	MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
	Xfr = gSum/gtPsi;

	// cout << "gSum: " << gSum << std::endl;
	
	
	// SBM mass matrix	
	HypreParMatrix Mmatp;
 	GridFunctionCoefficient cPs(&psi) ;	
	
	// use unique pointer to form ParBilinearForm
	std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(&fespace));
 	Mt->AddDomainIntegrator(new MassIntegrator(cPs)); 	
 	Mt->Assemble();
 	Mt->FormSystemMatrix(boundary_dofs, Mmatp);

	// Mmatp.Print("Mmatp_Original.txt", 0);

	
	HypreSmoother Mp_prec;
	CGSolver Mp_solver(MPI_COMM_WORLD);	

	Mp_solver.iterative_mode = false;
	Mp_solver.SetRelTol(1e-7);
	Mp_solver.SetAbsTol(0);
	Mp_solver.SetMaxIter(102);
	Mp_solver.SetPrintLevel(0);
	Mp_prec.SetType(HypreSmoother::Jacobi);
	Mp_solver.SetPreconditioner(Mp_prec);
	Mp_solver.SetOperator(Mmatp);

	HypreParMatrix *Tmatp;

	// SBM stiffness matrix
	ParGridFunction Dp(&fespace);		

	HypreParMatrix Kmatp;
	
// 	// force vector
	ParGridFunction Rxc(&fespace);	
	ParLinearForm Fct(&fespace);
	HypreParVector Fcb(&fespace);	

 	// create a Vector for CnP
	HypreParVector CpV0(&fespace), CpVn(&fespace), RHCp(&fespace);		
	
	int nDof = CpV0.Size();	
	
	// Vector of psi
	HypreParVector PsVc(&fespace);
	psi.GetTrueDofs(PsVc);		
	
	
// 	// 	// ============================
// 	// 	//    _____       ______ 
// 	// 	//   / ____|     |  ____|
// 	// 	//  | |     _ __ | |__   
// 	// 	//  | |    | '_ \|  __|  
// 	// 	//  | |____| | | | |____ 
// 	// 	//   \_____|_| |_|______|
// 	// 	// ============================
	
	ParGridFunction CnE(&fespace);
	double Ce0 = 0.001;						// initial value
	CnE = Ce0;	

	// SBM mass matrix	
	HypreParMatrix Mmate;	
	GridFunctionCoefficient cPe(&pse) ;	
	
	std::unique_ptr<ParBilinearForm> Me(new ParBilinearForm(&fespace)); 	
 	Me->AddDomainIntegrator(new MassIntegrator(cPe)); 	
 	Me->Assemble();
 	Me->FormSystemMatrix(boundary_dofs, Mmate);
	
	// Mmate.Print("Mmate_Original.txt", 0);

 	
	HypreSmoother Me_prec;
	CGSolver Me_solver(MPI_COMM_WORLD);	
	Me_solver.iterative_mode = false;
	Me_solver.SetRelTol(1e-7);
	Me_solver.SetAbsTol(0);
	Me_solver.SetMaxIter(100);
	Me_solver.SetPrintLevel(0);
	Me_prec.SetType(HypreSmoother::Jacobi);
	Me_solver.SetPreconditioner(Me_prec);
 	
	HypreParMatrix *TmatL, *TmatR;			// matrices for CN scheme
	
	
	// stiffness matrix
	ParGridFunction De(&fespace);				
	HypreParMatrix Kmate;

// 	// force vector
	ParGridFunction Rxe(&fespace);
	
// 	// defined for later
	ParLinearForm Fet(&fespace);		
	HypreParVector Feb(&fespace);	
	
	// // for imposing Neumann BC
	ParGridFunction PeR(&fespace);
	PeR = pse;
	PeR.Neg();
	GridFunctionCoefficient matCoef_R(&PeR) ;

	// PeR.Print(std::cout);

	// CnE.Print(std::cout);


	
 	// Vectors for CnE
	HypreParVector CeV0(&fespace), CeVn(&fespace), RHSe(&fespace);	
	
	// parameters used in the calculations
	double eCrnt = 0.0;
	double geCrnt = 0.0;
	double infx = 0.0;
	
	ParGridFunction CeT(&fespace);
	double CeC = 0.0;
	double CeAvg = 0.0;
	double gCeC = 0.0;			


// 	// ========================================
// 	// //               _     _____  
// 	// //              | |   |  __ \ 		//
// 	// //   _ __   ___ | |_  | |__) |
// 	// //  | '_ \ / _ \| __| |  ___/ 
// 	// //  | |_) | (_) | |_  | |     
// 	// //  | .__/ \___/ \__| |_|     
// 	// //  | |                       
// 	// //  |_|                       
// 	// ========================================

	ParGridFunction phP(&fespace);		// electropotential in particle
	ParGridFunction kap(&fespace);		// conductivity in particle
// 	ParGridFunction RpP(&fespace);		// reaction 
// 	ParGridFunction pP0(&fespace);		// values before iteration

	double BvP = 2.9395;
// 	BvP = 1.0;
	phP = BvP;	

	// stiffness matrix
	HypreParMatrix KmP;	

	CGSolver cgPP(MPI_COMM_WORLD);
	cgPP.SetRelTol(1e-7);
	cgPP.SetMaxIter(82);

// 	// force Vector
// 	ParLinearForm Fpt(&fespace);	
// 	HypreParVector Fpb(&fespace);

// 	HypreParVector Xs0(&fespace);


// 	// ========================================
// 	// // 	              _     ______ 
// 	// // 	             | |   |  ____|
// 	// // 	  _ __   ___ | |_  | |__   
// 	// // 	 | '_ \ / _ \| __| |  __|  
// 	// // 	 | |_) | (_) | |_  | |____ 
// 	// // 	 | .__/ \___/ \__| |______|
// 	// // 	 | |                       
// 	// // 	 |_|                       
// 	// ========================================
	
	double tc1 =(2*t_minus-1.0)/(2*t_minus*(1.0-t_minus));
	double tc2 = 1.0/(2*t_minus*(1.0-t_minus))*Cst1;
	double dffe;

	ParGridFunction phE(&fespace);		// electropot in electrolyte
	ParGridFunction Dmp(&fespace);		// D_minus_plus
	ParGridFunction kpl(&fespace);		// electrolyte conductivity
// 	ParGridFunction RpE(&fespace);		// reaction rate for electrolyte
// 	ParGridFunction pE0(&fespace);		// values before iteration

	double BvE = -1.0;
	phE = BvE;

	// stiffness matrix
	HypreParMatrix Kml;

	CGSolver cgPE(MPI_COMM_WORLD);
	cgPE.SetRelTol(1e-7);
	cgPE.SetMaxIter(80);
	 	
// 	// force vector
// 	ParLinearForm Flt(&fespace);
// 	HypreParVector Flb(&fespace);	
	
	// Laplace matrix
	HypreParMatrix Kdm;

	HypreParVector LpCe(&fespace), Xe0(&fespace);
// 	HypreParVector RHSl(&fespace);
	
// 	double Vcell = BvP - BvE;

	// reaction term
	ParGridFunction Rxn(&fespace);
	Rxn = 0.0;
	Rxn = AvP; Rxn *= 1.0e-8;

	// Rxn.Print(std::cout);	
	 
// 	// rate constants
// 	ParGridFunction dPHE(&fespace);		// voltage drop
	ParGridFunction Kfw(&fespace);		// forward reaction constant
	ParGridFunction Kbw(&fespace);		// backward rection constant
	ParGridFunction OCV(&fespace);		// open circuit voltage
	ParGridFunction i0C(&fespace);		// exchange current density	

	OCV = 0.0;
	i0C = 0.0;
	Kfw = 0.0;
	Kbw = 0.0;	

// 	double errP = 1.0;
// 	double errE = 1.0;
// 	int inlp = 0;
// 	double gErrP, gErrE;
// 	double dV, sgn;

// 	// containers used later.
	HypreParVector X1v(&fespace), B1v(&fespace);
	ParLinearForm B1t(&fespace);

// 	int cnt = 0;
// 	string stri;
// 	string fname;
		
// 	// 	===================================================================
// 	// 	//   _   _                      _                   _             
// 	// 	//  | | (_)                    | |                 (_)            
// 	// 	//  | |_ _ _ __ ___   ___   ___| |_ ___ _ __  _ __  _ _ __   __ _ 
// 	// 	//  | __| | '_ ` _ \ / _ \ / __| __/ _ \ '_ \| '_ \| | '_ \ / _` |
// 	// 	//  | |_| | | | | | |  __/ \__ \ ||  __/ |_) | |_) | | | | | (_| |
// 	// 	//   \__|_|_| |_| |_|\___| |___/\__\___| .__/| .__/|_|_| |_|\__, |
// 	// 	//                                     | |   | |             __/ |
// 	// 	//                                     |_|   |_|            |___/ 
// 	//  ===================================================================   
		
// 	int t = 0;
for (int t = 0; t < 10 + 1; t++){
// 	while ( Vcell > Vcut){
// //	while ( t<1 ){
	
	
// 		// // 	=============================
// 		// // 	   _____       _____  
// 		// // 	  / ____|     |  __ \ 		//
// 		// // 	 | |     _ __ | |__) |
// 		// // 	 | |    | '_ \|  ___/ 
// 		// // 	 | |____| | | | |     
// 		// // 	  \_____|_| |_|_|     
// 		// // 	==============================
	
		// Rxn.Print(std::cout);	
		Rxc = Rxn;
		Rxc /= rho;	
		GridFunctionCoefficient cAp(&Rxc) ;	

		// Rxc.Print(std::cout);
		
		// force term
		std::unique_ptr<ParLinearForm> Bc2(new ParLinearForm(&fespace));	
		// Use update () instead 
		Bc2->AddDomainIntegrator(new DomainLFIntegrator(cAp));	
		Bc2->Assemble();		
		// Move the contents of Bc2 into Fct
		Fct = std::move(*Bc2);
    
		// cout << "Bc2 in Original" << endl;
		// Bc2->Print(std::cout);		


		// diffusivity in the particles
		for (int vi = 0; vi < nV; vi++ ){
			// appendix equation A-19
			Dp(vi) = psi(vi)*(0.0277-0.084*CnP(vi) + 0.1003*CnP(vi)*CnP(vi))*1.0e-8;
			if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}

			// std::cout << "Dp[" << vi << "] = " << Dp(vi) << std::endl;
		}	
		GridFunctionCoefficient cDp(&Dp) ;	

		// std::cout << "cDp values in TimeStepCnP:" << std::endl;
		// 	for (int vi = 0; vi < fespace.GetTrueVSize(); ++vi) {
		// 		std::cout << (*cDp.GetGridFunction())(vi) << " ";
		// 	}
		// std::cout << std::endl;		
		
		// K matrix		
		std::unique_ptr<ParBilinearForm> Kc2(new ParBilinearForm(&fespace)); 
   		Kc2->AddDomainIntegrator(new DiffusionIntegrator(cDp));
   		Kc2->Assemble();
   		Kc2->FormLinearSystem(boundary_dofs, CnP, Fct, Kmatp, X1v, Fcb);
   		Fcb *= dt;

		// // Get the local data of the HypreParVector
		// double *Fxb_data = Fcb.GetData();

		// // Print each value of the vector
		// int size = Fcb.Size();
		// std::cout << "Fcb values in original:" << std::endl;
		// for (int i = 0; i < size; i++) {
		// 	std::cout << Fxb_data[i] << " ";
		// }
		// std::cout << std::endl;

		// T matrix
		Tmatp = Add(1.0, Mmatp, -dt, Kmatp);
		
		// vector of CnP				
		CnP.GetTrueDofs(CpV0);
		
		Tmatp->Mult(CpV0, RHCp);
		RHCp += Fcb;

		// // Get the local data of the HypreParVector
		// double *RHCp_data = RHCp.GetData();

		// // Print each value of the vector
		// int size = Fcb.Size();
		// std::cout << "Fcb values in original:" << std::endl;
		// for (int i = 0; i < size; i++) {
		// 	std::cout << RHCp_data[i] << " ";
		// }
		// std::cout << std::endl;		

		// time stepping
		Mp_solver.Mult(RHCp, CpVn) ;		

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

		// std::cout << "Updated CnP values:" << std::endl;
    	// CnP.Print(std::cout);
		
// 		delete Tmatp;


// 		// // ============================
// 		// //    _____       ______ 
// 		// //   / ____|     |  ____|
// 		// //  | |     _ __ | |__   
// 		// //  | |    | '_ \|  __|  
// 		// //  | |____| | | | |____ 
// 		// //   \_____|_| |_|______|
// 		// // ============================

		Rxe = Rxn;
		Rxe *= (-1.0*t_minus);	
		GridFunctionCoefficient cAe(&Rxe) ;

		Rxe.Print(std::cout);

		
		// total reaction
		eCrnt = 0.0;
		for (int ei = 0; ei < nE; ei++){	  
			Rxe.GetNodalValues(ei,VtxVal) ;
			val = 0.0;
			for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
			EAvg(ei) = val/nC;		 		  
			eCrnt += EAvg(ei)*EVol(ei);	
		} 	
		MPI_Allreduce(&eCrnt, &geCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		infx = geCrnt/L_w;

		// cout << "infx: " << infx << std::endl;
		
		
		// Neumann BC
		ConstantCoefficient nbcCoef(infx);

		ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);	

 		// force term
		std::unique_ptr<ParLinearForm> Be2(new ParLinearForm(&fespace));
		Be2->AddDomainIntegrator(new DomainLFIntegrator(cAe));
		Be2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m_nbcCoef), nbc_w_bdr);
		Be2->Assemble();
		// Move the contents of Be2 into Fet
		Fet = std::move(*Be2);

		// Fet.Print(std::cout);

				
		// salt diffusivity in the electrolyte					
		for (int vi = 0; vi < nV; vi++){
			// appendix equation A-21
			De(vi) = pse(vi)*D0*exp(-7.02-830*CnE(vi)+50000*CnE(vi)*CnE(vi));
			// std::cout << "De[" << vi << "] = " << De(vi) << std::endl;

			// std::cout << "De Diffusivity values in Original:" << std::endl;
			// for (int vi = 0; vi < std::min(nV, 10); ++vi) {
			// 	std::cout << De(vi) << " ";
			// }
    		// std::cout << std::endl;


		}
		GridFunctionCoefficient cDe(&De) ;

		// std::cout << "cDe values in Original CnE:" << std::endl;
		// 	for (int vi = 0; vi < fespace.GetTrueVSize(); ++vi) {
		// 		std::cout << (*cDe.GetGridFunction())(vi) << " ";
		// 	}
		// std::cout << std::endl;	
		
 		// K matrix		
		std::unique_ptr<ParBilinearForm> Ke2(new ParBilinearForm(&fespace)); 
   		Ke2->AddDomainIntegrator(new DiffusionIntegrator(cDe));
   		Ke2->Assemble();
   		Ke2->FormLinearSystem(boundary_dofs, CnE, Fet, Kmate, X1v, Feb);
   		Feb *= dt;	

		// // Get the local data of the HypreParVector
		// double *Feb_data = Feb.GetData();

		// // Print each value of the vector
		// int size1 = Feb.Size();
		// std::cout << "Feb values in original:" << std::endl;
		// for (int i = 0; i < size1; i++) {
		// 	std::cout << Feb_data[i] << " ";
		// }
		// std::cout << std::endl;
		
		// Crank-Nicolson matrices
		TmatR = Add(1.0, Mmate, -0.5*dt, Kmate);		
		TmatL = Add(1.0, Mmate,  0.5*dt, Kmate);		
		
		// vector of CnE				
		CnE.GetTrueDofs(CeV0);		
				
    	TmatR->Mult(CeV0, RHSe);
    	RHSe += Feb;
    	
    	// solver
		Me_solver.SetOperator(*TmatL);    	
    	
    	// time stepping
		Me_solver.Mult(RHSe, CeVn) ;
		
		// recover
		CnE.Distribute(CeVn);    	

		// check conservation of salt
		// if (t%500 == 0 && t > 0){
			CeC = 0.0;
			CeT = CnE;
			CeT *= pse;
			for (int ei = 0; ei < nE; ei++){
				CeT.GetNodalValues(ei,VtxVal) ;
				val = 0.0;
				for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
				EAvg(ei) = val/nC;	
				CeC += EAvg(ei)*EVol(ei) ;
			}
			MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
			// average CnE throughout electrolyte
			CeAvg = gCeC/gtPse;				
			
			// adjust CnE
			CnE -= (CeAvg-Ce0);
			MPI_Barrier(MPI_COMM_WORLD);
		// }	

		delete TmatR;
		delete TmatL;

		// std::cout << "Updated CnE values:" << std::endl;
    	// CnE.Print(std::cout);
		
		
// 		// ==============================================
// 		// // 	                      _   _             
// 		// // 	                     | | (_)            
// 		// // 	  _ __ ___  __ _  ___| |_ _  ___  _ __  
// 		// // 	 | '__/ _ \/ _` |/ __| __| |/ _ \| '_ \   //
// 		// // 	 | | |  __/ (_| | (__| |_| | (_) | | | |
// 		// // 	 |_|  \___|\__,_|\___|\__|_|\___/|_| |_|
// 		// ==============================================	


		// electrolyte conductivity and RHS	
		for (int vi = 0; vi < nV; vi++){
			dffe = exp(-7.02-830*CnE(vi)+50000*CnE(vi)*CnE(vi));
			Dmp(vi) = pse(vi)*tc1*D0*dffe;
			kpl(vi) = pse(vi)*tc2*D0*dffe*CnE(vi);

		}

		// std::cout << "Dmp = " << Dmp << std::endl;
		GridFunctionCoefficient cDm(&Dmp);

		// Laplace of CnE for the RHS
		std::unique_ptr<ParBilinearForm> Kl1(new ParBilinearForm(&fespace));
		Kl1->AddDomainIntegrator(new DiffusionIntegrator(cDm));
		Kl1->Assemble();
		Kl1->FormLinearSystem(boundary_dofs, phE, B1t, Kdm, X1v, B1v);		
	
		// Vector of CnE
		CnE.GetTrueDofs(CeVn) ;
		Kdm.Mult(CeVn, LpCe) ;

// 		// electrolyte conductivity and RHS		
		GridFunctionCoefficient cKe(&kpl) ;		
		std::unique_ptr<ParBilinearForm> Kl2(new ParBilinearForm(&fespace)); 		
		Kl2->AddDomainIntegrator(new DiffusionIntegrator(cKe));
		Kl2->Assemble();	

// 		// assign known values to the DBC nodes	
// 		ConstantCoefficient dbc_w_Coef(BvE);
		
// 		phE.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 		
		Kl2->FormLinearSystem(ess_tdof_list_w, phE, B1t, Kml, X1v, B1v);		
		
	// Solve the system using PCG with hypre's BoomerAMG preconditioner.
		//HypreBoomerAMG Mpe(Kml); //HypreBoomerAMG preconditioner causes memory leak issues
		//Mpe.SetPrintLevel(0);
		HypreSmoother Mpe;
		Mpe.SetType(HypreSmoother::Jacobi);
		cgPE.SetPreconditioner(Mpe);
		cgPE.SetOperator(Kml);
					
		
		// assign known values to the DBC nodes	
		// ConstantCoefficient dbc_e_Coef(BvP);			
		
		// particle conductivity
		// appendix equation A-20
		for (int vi = 0; vi < nV; vi++){
			kap(vi) = psi(vi)*(0.01929 + 0.7045*tanh(2.399*CnP(vi)) - \
				0.7238*tanh(2.412*CnP(vi)) - 4.2106e-6);
		}	

    	// std::cout << "kap = " << kap << std::endl;


		GridFunctionCoefficient cKp(&kap) ;
		
		// stiffness matrix for phP
		std::unique_ptr<ParBilinearForm> Kp2(new ParBilinearForm(&fespace));
		Kp2->AddDomainIntegrator(new DiffusionIntegrator(cKp));
		Kp2->Assemble();
		
		// project values to DBC nodes
		// phP.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 	
		Kp2->FormLinearSystem(ess_tdof_list_e, phP, B1t, KmP, X1v, B1v);			

		// Solve the system using PCG with hypre's BoomerAMG preconditioner.
		//HypreBoomerAMG Mpp(KmP);
		//Mpp.SetPrintLevel(0);
		HypreSmoother Mpp;
		Mpp.SetType(HypreSmoother::Jacobi);
		cgPP.SetPreconditioner(Mpp);
		cgPP.SetOperator(KmP);

		// std::cout << "AvB: " << AvB << std::endl;
		// std::cout << "CnP: " << CnP << std::endl;
		
		// rate constants and exchange current density at interface
		for (int vi = 0; vi < nV; vi++){
			if ( AvB(vi)*dh > 0.0 ){
				val = -0.2*(CnP(vi)-0.37)-1.559-0.9376*tanh(8.961*CnP(vi)-3.195); // check on this!
				i0C(vi) = pow(10.0,val)*1.0e-3;
				
				OCV(vi) = 1.095*CnP(vi)*CnP(vi) - 8.324e-7*exp(14.31*CnP(vi)) + \
					4.692*exp(-0.5389*CnP(vi));
					
				Kfw(vi) = i0C(vi)/(Frd*0.001  )*exp( alp*Cst1*OCV(vi)) ;	
				Kbw(vi) = i0C(vi)/(Frd*CnP(vi))*exp(-alp*Cst1*OCV(vi)) ;

				// std::cout << "Kbw[" << vi << "] = " << Kbw(vi) << std::endl;

			}
		}

		// std::cout << "Kbw = " << Kbw << std::endl;


// 		// convergence residuals
// 		gErrP = 1.0;
// 		gErrE = 1.0;
// 		inlp = 0;	

		// //  beginning of internal loop
		// while (gErrP > 1.0e-9 || gErrE > 1.0e-9 ){

		// 	// Butler-Volmer Equation for Reaction Rate
		// 	for (int vi = 0; vi < nV; vi++){
		// 		if ( AvB(vi)*dh > 0.0 ){
		// 			dPHE(vi) = phP(vi) - phE(vi);
		// 			Rxn(vi) = AvP(vi)*(Kfw(vi)*CnE(vi)*exp(-alp*Cst1*dPHE(vi)) - \
		// 			                   Kbw(vi)*CnP(vi)*exp( alp*Cst1*dPHE(vi)));
		// 		}
		// 	}


// 			// ========================================
// 			// //               _     _____  
// 			// //              | |   |  __ \ 		//
// 			// //   _ __   ___ | |_  | |__) |
// 			// //  | '_ \ / _ \| __| |  ___/ 
// 			// //  | |_) | (_) | |_  | |     
// 			// //  | .__/ \___/ \__| |_|     
// 			// //  | |                       
// 			// //  |_|                       
// 			// ========================================		
		
		
// 			// force vector
// 			RpP = Rxn;
// 			RpP *= Frd;	
// 			GridFunctionCoefficient cRp(&RpP);
		
// 			std::unique_ptr<ParLinearForm> Bp2(new ParLinearForm(&fespace));
// 			Bp2->AddDomainIntegrator(new DomainLFIntegrator(cRp));	
// 			Bp2->Assemble();		
// 			Fpt = std::move(*Bp2);		// Move the contents of Bp2 into Fpt

// 			// project values to DBC nodes
// 			phP.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 	
// 			Kp2->FormLinearSystem(ess_tdof_list_e, phP, Fpt, KmP, X1v, Fpb);			
			
// 			pP0 = phP;	
// 			pP0.GetTrueDofs(Xs0);
// 			cgPP.Mult(Fpb, Xs0);	
			
// 			// recover
// 			phP.Distribute(Xs0);   			

// 			for (int vi = 0; vi < nV; vi++){
// 				TmpF(vi) = pow(pP0(vi)-phP(vi),2)*psi(vi);
// 			}	
				
// 			errP = 0.0;
// 			for (int ei = 0; ei < nE; ei++){
// 				TmpF.GetNodalValues(ei,VtxVal) ;
// 				val = 0.0;
// 				for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
// 				EAvg(ei) = val/nC;						
// 				errP += EAvg(ei)*EVol(ei) ;
// 			}
// 			MPI_Allreduce(&errP, &gErrP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

// 			gErrP /= gtPsi;
// 			gErrP = pow(gErrP, 0.5);
			
			
// 			// ========================================
// 			// // 	              _     ______ 
// 			// // 	             | |   |  ____|
// 			// // 	  _ __   ___ | |_  | |__   
// 			// // 	 | '_ \ / _ \| __| |  __|  
// 			// // 	 | |_) | (_) | |_  | |____ 
// 			// // 	 | .__/ \___/ \__| |______|
// 			// // 	 | |                       
// 			// // 	 |_|                       
// 			// ========================================
		
		
// 			RpE = Rxn;
// 			RpE.Neg();	
// 			GridFunctionCoefficient cRe(&RpE) ;				
			
// 			std::unique_ptr<ParLinearForm> Bl2(new ParLinearForm(&fespace));
// 			Bl2->AddDomainIntegrator(new DomainLFIntegrator(cRe));	
// 			Bl2->Assemble();		
// 			Flt = std::move(*Bl2);		// Move the contents of Bl2 into Flt
		
// 			phE.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 		
// 			Kl2->FormLinearSystem(ess_tdof_list_w, phE, Flt, Kml, X1v, Flb);		
						
// 			RHSl = Flb;
// 			RHSl += LpCe;

// 			pE0 = phE;
// 			pE0.GetTrueDofs(Xe0);
// 			cgPE.Mult(RHSl, Xe0);   
			
// 			// recover
// 			phE.Distribute(Xe0);		

// 			for (int vi = 0; vi < nV; vi++){
// 				TmpF(vi) = pow(pE0(vi)-phE(vi),2)*pse(vi);
// 			}	

// 			errE = 0.0;	
// 			for (int ei = 0; ei < nE; ei++){
// 				TmpF.GetNodalValues(ei,VtxVal) ;
// 				val = 0.0;
// 				for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
// 				EAvg(ei) = val/nC;	
// 				errE += EAvg(ei)*EVol(ei) ;					
// 			}	
// 			MPI_Allreduce(&errE, &gErrE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
// 			gErrE /= gtPse;
// 			gErrE = pow(gErrE,0.5);
			
// 			inlp += 1;

// 		} // end of internal loop


// 		// total reaction current
// 		sCrnt = 0.0;
// 		for (int ei = 0; ei < nE; ei++){
// 			Rxn.GetNodalValues(ei,VtxVal) ;
// 			val = 0.0;
// 			for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
// 			EAvg(ei) = val/nC;
// 			sCrnt += EAvg(ei)*EVol(ei) ;
// 		} 			
// 		MPI_Allreduce(&sCrnt, &gCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
					
// 		// adjust BvP for constant C rate loading
// 		sgn = copysign(1, gTrgI - gCrnt);
// 		dV = dt*Vsr*sgn;
// 		BvP -= dV;
// 		phP -= dV;
		
// 		Vcell = BvP - BvE;

// // // 		cout << myid << " " << sgn << " " << gTrgI-gCrnt << endl;
		
// 		tm = tm + dt;

		// t += 1;

// 		if ( Xfr > 0.3+cnt*0.005 ){
// 			// output values to a text file
// 			if (myid == 1){
// 				ofstream myfile;
// 				myfile.open ("Output_2D_D01_01.txt", ios::app);
// 				myfile << cnt << "  " << t << "  " << tm << "  " << Xfr << "  " << Vcell << "  " \
// 					<< "  " << gCrnt << "  " << gTrgI << "  " << CeAvg << endl;
// 				myfile.close();
// 			}
			
// 			stri = to_string(1000+cnt);
		
// 			fname = "Output/CnP_80x90_D01_" + stri;			
// 			const char *test1 = fname.c_str();
// 			CnP.Save(test1);

// 			fname = "Output/CnE_80x90_D01_" + stri;			
// 			const char *test2 = fname.c_str();
// 			CnE.Save(test2);
			
// 			fname = "Output/phP_80x90_D01_" + stri;			
// 			const char *test3 = fname.c_str();
// 			phP.Save(test3);			

// 			fname = "Output/phE_80x90_D01_" + stri;			
// 			const char *test4 = fname.c_str();
// 			phE.Save(test4);
						
// 			cnt += 1;
			
// 			MPI_Barrier(MPI_COMM_WORLD);
// 		}
		
// 		if ( Vcell < Vcut ){
// 			// output values to a text file
// 			if (myid == 1){
// 				ofstream myfile;
// 				myfile.open ("Output_2D_D01_01.txt", ios::app);
// 				myfile << cnt << "  " << t << "  " << tm << "  " << Xfr << "  " << Vcell << "  " \
// 					<< "  " << gCrnt << "  " << gTrgI << "  " << CeAvg << endl;
// 				myfile.close();
// 			}
			
// 			stri = to_string(1000+cnt);
		
// 			fname = "Output/CnP_80x90_D01_" + stri;			
// 			const char *test1 = fname.c_str();
// 			CnP.Save(test1);

// 			fname = "Output/CnE_80x90_D01_" + stri;			
// 			const char *test2 = fname.c_str();
// 			CnE.Save(test2);
			
// 			fname = "Output/phP_80x90_D01_" + stri;			
// 			const char *test3 = fname.c_str();
// 			phP.Save(test3);			

// 			fname = "Output/phE_80x90_D01_" + stri;			
// 			const char *test4 = fname.c_str();
// 			phE.Save(test4);
						
// 			cnt += 1;
						
// 		}
		
		
// 	if (t%200 == 1 && myid == 1){cout << t << " - " << Xfr << " - " << tm << " - " << \
// 		Vcell << endl;}
// // 	if (myid == 1 ){cout << t << "  " << Xfr << "  " << tm << "  " << Vcell << endl;}
		
// 	} // time iteration loop


// // 	fname = "Output/CnP_80x90_P09_" + stri;			
// // 	const char *test1 = fname.c_str();
// // 	CnP.Save(test1);	
// // 	
// // 	fname = "Output/CnE_80x90_P09_" + stri;			
// // 	const char *test2 = fname.c_str();
// // 	CnE.Save(test2);	
// // 	
// // 	fname = "Output/phP_80x90_P09_" + stri;			
// // 	const char *test3 = fname.c_str();
// // 	phP.Save(test3);			
// // 
// // 	fname = "Output/phE_80x90_P09_" + stri;			
// // 	const char *test4 = fname.c_str();
// // 	phE.Save(test4);	
	
	
	
//    return 0;
	}

	// int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the MPI rank

    // std::string file_name = "CnP_solution_original." + std::to_string(rank) + ".gf";  // MFEM's GridFunction format (.gf)
    
    // std::ofstream ofs(file_name.c_str());
    // if (ofs.is_open()) {
    //     CnP.Save(ofs);  // Use '->' because CnE is a pointer
    //     ofs.close();
    // } else {
    //     mfem::mfem_error("Error opening file to save CnE.");
    // }

	// int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the MPI rank

    // std::string file_name = "CnE_solution_original." + std::to_string(rank) + ".gf";  // MFEM's GridFunction format (.gf)
    
    // std::ofstream ofs(file_name.c_str());
    // if (ofs.is_open()) {
    //     CnE.Save(ofs);  // Use '->' because CnE is a pointer
    //     ofs.close();
    // } else {
    //     mfem::mfem_error("Error opening file to save CnE.");
    // }
}

