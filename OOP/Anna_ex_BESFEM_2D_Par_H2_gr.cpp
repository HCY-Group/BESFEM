#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "mpi.h"
#include <chrono>


#include <cmath>

// #include <stdio.h>
	// const char *mesh_file = "../Inputs/Mesh_96x96_P01.mesh";
	// const char *dsF_file =  "../Inputs/dsF_96x96_P01.txt";

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{

    // Start measuring the program execution time
    using namespace std::chrono;
    auto program_start = high_resolution_clock::now();

    // Initialize MPI for parallel processing and HYPRE for solver setup
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();
 
{  
   int myid = Mpi::WorldRank();
   
	// 1. Parse command line options.
	
	// const char *mesh_file = "../Inputs/Mesh_81x81x6_disk.mesh";
	// // const char *dsF_file =  "../Inputs/DSTFB_81x81x6.txt";

	// const char *mesh_file = "../Inputs/HCY_Mesh_80x80x16_01.mesh";
	// const char *dsF_file =  "../Inputs/HCY_dsF_81x81x17_01.txt";

	// const char *mesh_file = "../Inputs/Mesh_80x80x5_3D_disk.mesh";
	// const char *dsF_file =  "../Inputs/dsF_80x80x5_3D_disk.txt";

	// const char *mesh_file = "../Inputs/HCY_Mesh_10x10x2_01.mesh";
	// const char *dsF_file =  "../Inputs/HCY_dsF_11x11x3_01.txt";

	const char *mesh_file = "../Inputs/disk_Mesh_80x80x6.mesh";
	const char *dsF_file =  "../Inputs/disk_dsF_81x81x7.txt";
	


// 	bool visualization = true;

	int AMR = 0;
	bool ReMsh = false; 	

	int order = 1;

	double dh = 3.25e-5;
	// double dh = 0.4e-4;
	// double dh = 2.6e-4;

	double zeta = 1.0;				    // interfacial thickness
	double thres = 1.0e-2;					// AvP musk
	double eps = 1.0e-6;					// var-epsilon			
	double dt = 0.0105625 * 1.0;	// time step
	// double dt = 0.00792188;
	// double dt = 0.00528125;
	double tm = 0.0;		
	double gc = 3.3800e-10 *2;			// gradient coefficient				
	
	double t_minus = 7.619047619047619e-01; // transference number
	double D0 = 0.00489; 					// base diffusivity	
	double Frd = 96485.3365; 				// Faraday constant
	double Cst1 = 1.6021766e-19/(1.3806488e-23*300.0);
	double alp = 0.5;	
	
	double rho = 0.0312;		 	// Li site density	
	double Cr = 0.5;				// C-rate
	// double Cr = 0.25; 			// C-rate HCY 
	double Vsr = 0.009466;			// voltage scanning rate
	double Vcut = 0.0; 				// cut-off voltage
	


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

	// cout << "Mesh dimension = " << gmesh.Dimension() << endl;
	// cout << "Boundary attributes in mesh: " << gmesh.bdr_attributes << endl;
	

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

	double min_dsF = gDsF.Min();
	double max_dsF = gDsF.Max();
	// cout << "dsF range: " << min_dsF << " to " << max_dsF << endl;

	
	// AMR =========================================
	int nCg = pow(2,gmesh.Dimension());
	Array<double> VertVal(nCg) ;
	double vSum = 0.0;
	
	// levels of refinement
	if (AMR > 0){
		ReMsh = true;
		
		// create an array to store vertex labels of each element
		double Ds2Br ;		// distance of an element to the internal boundary	
		
		int Lv = AMR; 
		Array<double> LvFc(Lv) ; // criteria for refinement
		// for (int i = 0; i < Lv; i++){
		// 	LvFc[i] = dh * pow(2.0, -i);  

		// }

		LvFc[0] = -1.5;
		LvFc[1] =  -0.8;
		LvFc[2] =  -0.5;
		// LvFc[3] =  -2.60;
		// LvFc = [-0.8, -0.4, 0.0]
	
		// create an array to store the labels of elements needed to be refined.
		Array<int> refinement_list;
	
		for (int Lv = 0; Lv < AMR; Lv++){
			
			int refine_count = 0;

			for (int e = 0; e < gmesh.GetNE(); e++){
				// get the distance values at the 4 vertices
				gDsF.GetNodalValues(e,VertVal) ;
		
				// distance from element center to the internal boundary
				vSum = 0.0;
				for (int vt = 0; vt < nCg; vt++){
					vSum += VertVal[vt];
				}
				Ds2Br = vSum/nCg;	
		
				if (Ds2Br > LvFc[Lv])
				{
					refinement_list.Append(e);
					refine_count++;
				}
			}

			cout << "Refinement level " << Lv << ": refining " << refine_count << " elements." << endl;


			// mesh refinement
			gmesh.GeneralRefinement(refinement_list);
	
			// update relevant objects 
			gFespace.Update();
			gDsF.Update() ;
		
			// clear refinement list
			refinement_list.DeleteAll() ;
		}	
	}
	// =========================================

	// H1_FECollection gfec(order, gmesh.Dimension());
	// FiniteElementSpace gfespace(&gmesh, &gfec);

	// GridFunction gDsF(&gfespace);	// global distance function
	// gDsF.ProjectGridFunction(gDsF_in); // Project the global distance function

	// west boundary size
	Vector Rmin, Rmax;
	gmesh.GetBoundingBox(Rmin,Rmax);
	// double L_w = Rmax(1) - Rmin(1);
	// double L_w = (Rmax(1) - Rmin(1)) + 2*(Rmax(0) - Rmin(0));
	// double L_w = (Rmax(1) - Rmin(1))*(Rmax(2) - Rmin(2));
	double L_w = (Rmax(1) - Rmin(1))*(Rmax(2) - Rmin(2)); //3D
	
	// local (parallel) mesh 
	ParMesh pmesh(MPI_COMM_WORLD, gmesh);
	pmesh.Save("Results/Pmesh_96x96_P01");

	
	int nV = pmesh.GetNV();					// number of vertices
	int nE = pmesh.GetNE();					// number of elements
	int nC = pow(2,pmesh.Dimension());		// number of corner vertices
	
	Array<double> VtxVal(nC) ;				// values of corner vertices
	
	// element area
	Vector EVol(nE);
	for (int ei = 0; ei < nE; ei++){
		EVol(ei) = pmesh.GetElementVolume(ei);	
	} 	
	
	
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
		gmesh.GetElementVertices(gei,gVTX);
		pmesh.GetElementVertices(ei,VTX);
	
		for (int vi = 0; vi < nC; vi++){
			dsF(VTX[vi]) = gDsF(gVTX[vi]); // Convert distance
		}			
	}	
	
	// Local domain parameters
	ParGridFunction psi(&fespace);
	ParGridFunction pse(&fespace);	
	ParGridFunction AvP(&fespace);	

	// interpolate domain parameter from distance function
    for(int vi = 0; vi < nV; vi++){
        
		// double dsF_um = dsF(vi);
		// double dsF_cm = dsF_um * 1.0e-4;  // μm → cm
		// double arg = dsF_cm / (zeta * dh);

		// psi(vi) = 0.5 * (1.0 + tanh(arg));
		// pse(vi) = 1.0 - psi(vi);
		// AvP(vi) = -(pow(tanh(arg),2) - 1.0) / (2.0 * zeta * dh);

		// dsF(vi) *= 1.0e-4;  // Convert distance from μm to cm
		// psi(vi) = 0.5*(1.0 + tanh(dsF(vi)/(zeta*dh)));  
        // pse(vi) = 1.0 - psi(vi);
        // AvP(vi) = -(pow(tanh(dsF(vi)/(zeta*dh)),2) - 1.0)/(2*zeta*dh);

		psi(vi) = 0.5*(1.0 + tanh(dsF(vi)/(zeta)));  
        pse(vi) = 1.0 - psi(vi);
        AvP(vi) = -(pow(tanh(dsF(vi)/(zeta)),2) - 1.0)/(2*zeta);
        
        // if (psi(vi) < eps){psi(vi) = eps;}   
        // if (pse(vi) < eps){pse(vi) = eps;}  
		
		if (psi(vi) < 0){psi(vi) = 0;}   
		if (psi(vi) > 1){psi(vi) = 1;}  
        if (pse(vi) < 0){pse(vi) = 0;}  
		if (pse(vi) > 1){pse(vi) = 1;}       
    } 	

	// double dh_cm = 3.25e-4;         // 3.25 μm = 3.25e-4 cm
	// double dx_cm = 3.25e-4;         // mesh spacing in cm

	// for (int vi = 0; vi < nV; vi++)
	// {
	// 	dsF(vi) *= 1.0e-4; // Convert distance from μm to cm

	// 	double arg = dsF(vi) / (zeta * dh_cm);
	// 	double tanh_arg = tanh(arg);

	// 	// AvP in cm^-1
	// 	double AvP_cm_inv = (1.0 - tanh_arg * tanh_arg) / (2.0 * zeta * dh_cm);

	// 	// Convert AvP to cm^-2 for use in reaction source term
	// 	AvP(vi) = AvP_cm_inv / dx_cm;

	// 	psi(vi) = 0.5 * (1.0 + tanh_arg);
	// 	pse(vi) = 1.0 - psi(vi);

	// 	if (psi(vi) < eps) { psi(vi) = eps; }
	// 	if (pse(vi) < eps) { pse(vi) = eps; }
	// }

	psi += 1.0e-6;
	pse += 1.0e-6;

	// output psi
	psi.Save("Results/psi_96x96_P01");
	

	ParGridFunction AvB(&fespace);
	AvB = AvP;
	for (int vi = 0; vi < nV; vi++){
		 // remove small values of AvP  
		// if (AvP(vi)*dh < 1.0e-3 ){AvP(vi) = 0.0;} 
		// if (AvB(vi)*dh < 1.0e-6){AvB(vi) = 0.0;} 
		if (AvP(vi) < 1.0e-2 ){AvP(vi) = 0.0;} 
		if (AvB(vi) < 1.0e-6){AvB(vi) = 0.0;} 
		 		 
	}
	AvP /= dh;

	// const double max_AvP = 13120.0;

	// for (int vi = 0; vi < AvP.Size(); ++vi)
	// {
	// 	if (AvP(vi) > max_AvP)
	// 	{
	// 		AvP(vi) = max_AvP;
	// 	}
	// }
 	AvP.Save("Results/AvP_96x96_P01");	
	AvB.Save("Results/AvB_96x96_P01");


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
	double trgI = tPsi*rho*(1.0-0.0)/(3600.0/Cr);
	double gTrgI = 0.0;
	MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	double sCrnt = 0.0;	
	double gCrnt = 0.0;	

	cout << "Total psi = " << gtPsi << endl;
	cout << "Total pse = " << gtPse << endl;
	cout << "Target current = " << gTrgI << endl;

	pmesh.PrintInfo();

	// // Some containers that will be used later.
	// Array<int> boundary_dofs;					// nature boundary
	
	// // Neumann BC on the west boundary. CnE
	// Array<int> nbc_w_bdr(pmesh.bdr_attributes.Max());
	// nbc_w_bdr = 0; nbc_w_bdr[1] = 1; nbc_w_bdr[2] = 1; nbc_w_bdr[3] = 1;

	// // Dirichlet BC on the east boundary. phP
	// Array<int> dbc_e_bdr(pmesh.bdr_attributes.Max());
	// dbc_e_bdr = 0; dbc_e_bdr[0] = 1;
	// // use dbc_e_bdr array to extract all node labels of Dirichlet BC
	// Array<int> ess_tdof_list_e(0);			
	// fespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);
	
	// // Dirichlet BC on the west boundary. phE
	// Array<int> dbc_w_bdr(pmesh.bdr_attributes.Max());
	// dbc_w_bdr = 0; dbc_w_bdr[1] = 1; dbc_w_bdr[2] = 1; dbc_w_bdr[3] = 1;
	// // use dbc_w_bdr array to extract all node labels of Dirichlet BC
	// Array<int> ess_tdof_list_w(0);			
	// fespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

	
	// Some containers that will be used later.
	Array<int> boundary_dofs;					// nature boundary
	
	// Neumann BC on the west boundary. CnE
	Array<int> nbc_w_bdr(pmesh.bdr_attributes.Max());
	nbc_w_bdr = 0;
	// nbc_w_bdr[4] = 1; // top boundary
	// nbc_w_bdr[5] = 1; // bottom boundary
	nbc_w_bdr[2] = 1;	

	// Dirichlet BC on the west boundary. phE 
	Array<int> dbc_w_bdr(pmesh.bdr_attributes.Max());
	dbc_w_bdr = 0; 
	// dbc_w_bdr[0] = 1;
	dbc_w_bdr[2] = 1; // east dirch for electrolyte
	// dbc_w_bdr[3] = 1;
	// dbc_w_bdr[1] = 1;

		// CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phE(0:ny+1,0:nx+1,0:nz+1) )
		// phE(0:ny+1,0,0:nz+1) = BvE
		// phE(0:ny+1,nx+1,0:nz+1) = phE(0:ny+1,nx,0:nz+1)
		// dirichlet BC on the west phE


	Array<int> ess_tdof_list_w(0);			
	fespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

	// Dirichlet BC on the east phP
	Array<int> dbc_e_bdr(pmesh.bdr_attributes.Max());
	dbc_e_bdr = 0; 
	dbc_e_bdr[0] = 1; // this is on the west dirch for particle
	// dbc_e_bdr[2] = 1; // this is on the east
	// dbc_e_bdr[3] = 1;
	// dbc_e_bdr[1] = 1;

		// CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phP(0:ny+1,0:nx+1,0:nz+1) )	
		// phP(0:ny+1,nx+1,0:nz+1) = BvP
		// phP(0:ny+1,0,0:nz+1) = phP(0:ny+1,1,0:nz+1)
		// dirichlet BC on the east phP	

	Array<int> ess_tdof_list_e(0);
	fespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

		// Get psi values at all true DOFs
	// Vector psi_vals;
	// psi.GetTrueDofs(psi_vals);

	// Get raw essential DOFs
	// Array<int> raw_ess_phP, raw_ess_phE;
	// fespace.GetEssentialTrueDofs(dbc_e_bdr, raw_ess_phP);
	// fespace.GetEssentialTrueDofs(dbc_w_bdr, raw_ess_phE);

	// // Filter based on psi
	// for (int i = 0; i < raw_ess_phP.Size(); i++) {
	// 	int dof = raw_ess_phP[i];
	// 	if (psi_vals[dof] > 0.5)
	// 		ess_tdof_list_e.Append(dof);  // graphite only
	// }

	// for (int i = 0; i < raw_ess_phE.Size(); i++) {
	// 	int dof = raw_ess_phE[i];
	// 	if (psi_vals[dof] < 0.5)
	// 		ess_tdof_list_w.Append(dof);  // electrolyte only
	// }

	// cout << "phP essential tdofs: " << ess_tdof_list_e.Size() << endl;
	// cout << "phE essential tdofs: " << ess_tdof_list_w.Size() << endl;


			
	// for (int i = 0; i < pmesh.GetNBE(); i++) {
    // int attr = pmesh.GetBdrAttribute(i);
    // cout << "Boundary element " << i << " has attribute " << attr << endl;
	// }
	
	// delete global mesh because it is longer used
	gmesh.Clear();
		
	MPI_Barrier(MPI_COMM_WORLD);
	
	// containers used later.
	HypreParVector X1v(&fespace), B1v(&fespace);
	ParLinearForm B1t(&fespace);
		
	// // 	=============================
	// // 	   _____       _____  
	// // 	  / ____|     |  __ \ 		//
	// // 	 | |     _ __ | |__) |
	// // 	 | |    | '_ \|  ___/ 
	// // 	 | |____| | | | |     
	// // 	  \_____|_| |_|_|     
	// // 	==============================

	// initial condition
	ParGridFunction CnP(&fespace);
	double Cp0 = 2.02e-2;					// initial value
	CnP = Cp0;	
	
	double Ctmp;
		
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
	
	// std::cout << "Initial degree of lithiation Xfr = " << Xfr << std::endl;
	
	// SBM mass matrix	
	HypreParMatrix Mmatp;
 	GridFunctionCoefficient cPs(&psi) ;	
	// ConstantCoefficient c1(1.0);

	
	// use unique pointer to form ParBilinearForm
	std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(&fespace));
	// Mt->AddDomainIntegrator(new MassIntegrator(c1));
 	Mt->AddDomainIntegrator(new MassIntegrator(cPs)); 	
 	Mt->Assemble();
 	Mt->FormSystemMatrix(boundary_dofs, Mmatp);
	
	HypreSmoother Mp_prec;
	Mp_prec.SetType(HypreSmoother::Jacobi);

	// HypreBoomerAMG Mp_prec;
	// Mp_prec.SetPrintLevel(0);
	
	CGSolver Mp_solver(MPI_COMM_WORLD);	
	// Mp_solver.iterative_mode = false;
	Mp_solver.SetRelTol(1e-6);
	Mp_solver.SetAbsTol(0);
	Mp_solver.SetMaxIter(102);
	Mp_solver.SetPrintLevel(0);
	Mp_solver.SetPreconditioner(Mp_prec);
	Mp_solver.SetOperator(Mmatp);

	// SBM stiffness matrix
	ParGridFunction Mob(&fespace);	
		
	// Mob = 1.0e-12; // avg 6.0e-12
	// Mob *= psi;
	
	GridFunctionCoefficient cDp(&Mob) ;	

	// GridFunctionCoefficient cMob(&Mob);
	// GridFunctionCoefficient cPsi(&psi);
	// ProductCoefficient mob_with_psi(cMob, cPsi);
	
	HypreParMatrix Kmatp;
	std::unique_ptr<ParBilinearForm> Kc2(new ParBilinearForm(&fespace)); 
	// Kc2->AddDomainIntegrator(new DiffusionIntegrator(mob_with_psi));
	// Kc2->AddDomainIntegrator(new DiffusionIntegrator(cMob));
	Kc2->AddDomainIntegrator(new DiffusionIntegrator(cDp));
	Kc2->Assemble();
	
	// gradient energy term, Lap, Mat
	// ConstantCoefficient varE(0.32 * (dh * dh) * 2); // eps in fortran interfacial thickness
	// try varE(6.7600e-8) new and varE(6.76e-10) original; 
	// ConstantCoefficient varE(6.7600e-8); // eps in fortran
	ConstantCoefficient  varE(gc/pow(dh,pmesh.Dimension())); // <== note here 
	ParGridFunction Mub(&fespace);

	// double varE_val = 1.2 * (dh * dh);
	// ConstantCoefficient varE(varE_val);
	// std::cout << "varE = " << varE_val << std::endl;
			
	HypreParMatrix KmX;
	std::unique_ptr<ParBilinearForm> Kx2(new ParBilinearForm(&fespace)); 
	Kx2->AddDomainIntegrator(new DiffusionIntegrator(varE));
	Kx2->Assemble();		
	
	// force vector
	ParGridFunction Rxc(&fespace);
	GridFunctionCoefficient cAp(&Rxc) ;		
	
	ParLinearForm Fct(&fespace);
	HypreParVector Fcb(&fespace);	

	std::unique_ptr<ParLinearForm> Bc2(new ParLinearForm(&fespace));		
	Bc2->AddDomainIntegrator(new DomainLFIntegrator(cAp));	
	Bc2->Assemble();		
	Fct = *Bc2;		// Move the contents of Bc2 into Fct	
	
	// varE nabla^2 phi	
	Kx2->FormLinearSystem(boundary_dofs, CnP, Fct, KmX, X1v, Fcb);	
	
	// KmX.Print("Kmx.txt", 0);
	
	// nabla cdot Mob nabla Mub	
	Kc2->FormLinearSystem(boundary_dofs, Mub, Fct, Kmatp, X1v, Fcb);
	
	
 	// create a Vector for CnP
	HypreParVector CpV0(&fespace), CpVn(&fespace), RHCp(&fespace);		
	
	int nDof = CpV0.Size();	

	HypreParVector Lp1(&fespace), Lp2(&fespace), MuV(&fespace);	
	
	// Vector of psi
	HypreParVector PsVc(&fespace);
	psi.GetTrueDofs(PsVc);		

	
	// tables 
	int pnt1 = 101;
	vector<double> Ticks(pnt1);
	ifstream tck;
	std::string dX_File = "../Inputs/C_Li_X_101.txt";
	tck.open(dX_File);
	for (int i = 0; i < pnt1; i++){
		tck >> Ticks[i];		
	}
	tck.close();
	double d_X = Ticks[2] - Ticks[1];		
	
	// chemical potential
	vector<double> chmPot(pnt1);
	ifstream dF_val;
	std::string dF_File = "../Inputs/C_Li_M6_101.txt";
	dF_val.open(dF_File);
	for (int i = 0; i < pnt1; i++){
		dF_val >> chmPot[i];	
	}
	dF_val.close();
	
	// mobility
	vector<double> Mobil(pnt1);
	ifstream dM_val;
	std::string dM_File = "../Inputs/C_Li_Mb5_101.txt";
	dM_val.open(dM_File);
	for (int i = 0; i < pnt1; i++){
		dM_val >> Mobil[i];
		Mobil[i] *= 100.0 * 2.0 / 3.0;  // Apply Fortran scaling

	}
	dM_val.close();

	// for (int i = 0; i < 101; ++i)
    // std::cout << "Mobil[" << i << "] = " << Mobil[i] << std::endl;

	
	// OCV
	vector<double> OCVtb(pnt1);
	ifstream dO_val;
	std::string dO_File = "../Inputs/C_Li_O3_101.txt";
	dO_val.open(dO_File);
	for (int i = 0; i < pnt1; i++){
		dO_val >> OCVtb[i];
	}
	dO_val.close();
	
	// i0
	vector<double> i0tb(pnt1);
	ifstream dJ_val;
	std::string dJ_File = "../Inputs/C_Li_J2_101.txt";
	dJ_val.open(dJ_File);
	for (int i = 0; i < pnt1; i++){
		dJ_val >> i0tb[i];
	}
	dJ_val.close();	
	
	
	// 	// ============================
	// 	//    _____       ______ 
	// 	//   / ____|     |  ____|
	// 	//  | |     _ __ | |__   
	// 	//  | |    | '_ \|  __|  
	// 	//  | |____| | | | |____ 
	// 	//   \_____|_| |_|______|
	// 	// ============================
	
	ParGridFunction CnE(&fespace);
	double Ce0 = 0.001005;						// initial value
	CnE = Ce0;	

	// SBM mass matrix	
	HypreParMatrix Mmate;	
	GridFunctionCoefficient cPe(&pse) ;	
	
	std::unique_ptr<ParBilinearForm> Me(new ParBilinearForm(&fespace)); 	
 	Me->AddDomainIntegrator(new MassIntegrator(cPe)); 	
 	Me->Assemble();
 	Me->FormSystemMatrix(boundary_dofs, Mmate);
 	
	HypreSmoother Me_prec;
	Me_prec.SetType(HypreSmoother::Jacobi);
	// HypreBoomerAMG Me_prec;
	// Me_prec.SetPrintLevel(0);
	
	CGSolver Me_solver(MPI_COMM_WORLD);	
	// Me_solver.iterative_mode = false;
	Me_solver.SetRelTol(1e-6);
	Me_solver.SetAbsTol(0);
	Me_solver.SetMaxIter(100);
	Me_solver.SetPrintLevel(0);
	Me_solver.SetPreconditioner(Me_prec);
 	
	HypreParMatrix *TmatL, *TmatR;			// matrices for CN scheme
	
	// stiffness matrix
	ParGridFunction De(&fespace);	
	GridFunctionCoefficient cDe(&De) ;
				
	HypreParMatrix Kmate;
	std::unique_ptr<ParBilinearForm> Ke2(new ParBilinearForm(&fespace)); 
	Ke2->AddDomainIntegrator(new DiffusionIntegrator(cDe));
	Ke2->Assemble();

	double infx = 0.0;
	// for imposing Neumann BC
	ParGridFunction PeR(&fespace);
	PeR = pse;
	PeR.Neg();	
	GridFunctionCoefficient matCoef_R(&PeR);
	
	ConstantCoefficient nbcCoef(infx);
	ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

	// force vector
	ParGridFunction Rxe(&fespace);
	GridFunctionCoefficient cAe(&Rxe) ;	
	
	// defined for later
	ParLinearForm Fet(&fespace);		
	HypreParVector Feb(&fespace);

	// Initialize linear form
	std::unique_ptr<ParLinearForm> Be2(new ParLinearForm(&fespace));
	Be2->AddDomainIntegrator(new DomainLFIntegrator(cAe));
	Be2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m_nbcCoef), nbc_w_bdr);
	Be2->Assemble();
		
	Fet = *Be2;	
	
	Ke2->FormLinearSystem(boundary_dofs, CnE, Fet, Kmate, X1v, Feb);
   		
 	// Vectors for CnE
	HypreParVector CeV0(&fespace), CeVn(&fespace), RHSe(&fespace);	
	
	// parameters used in the calculations
	double eCrnt = 0.0;
	double geCrnt = 0.0;

	ParGridFunction CeT(&fespace);
	double CeC = 0.0;
	double CeAvg = 0.0;
	double gCeC = 0.0;			


	// ========================================
	// //               _     _____  
	// //              | |   |  __ \ 		//
	// //   _ __   ___ | |_  | |__) |
	// //  | '_ \ / _ \| __| |  ___/ 
	// //  | |_) | (_) | |_  | |     
	// //  | .__/ \___/ \__| |_|     
	// //  | |                       
	// //  |_|                       
	// ========================================

	ParGridFunction phP(&fespace);		// electropotential in particle
	ParGridFunction pP0(&fespace);		// values before iteration
		
	ParGridFunction kap(&fespace);		// conductivity in particle
	kap = psi;
	kap *= 3.3;
	GridFunctionCoefficient cKp(&kap) ;
	
	ParGridFunction RpP(&fespace);		// reaction 
	GridFunctionCoefficient cRp(&RpP);

	double BvP = -0.1;
	phP = BvP;	

	// stiffness matrix
	HypreParMatrix KmP;	

	std::unique_ptr<ParBilinearForm> Kp2(new ParBilinearForm(&fespace));
	Kp2->AddDomainIntegrator(new DiffusionIntegrator(cKp));
	Kp2->Assemble();

	// // project values to DBC nodes
	// ConstantCoefficient dbc_e_Coef(BvP);	
	// phP.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 	
	// Kp2->FormLinearSystem(ess_tdof_list_e, phP, B1t, KmP, X1v, B1v);	

	// project values to DBC nodes
	ConstantCoefficient dbc_e_Coef(BvP);	
	phP.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 	
	Kp2->FormLinearSystem(ess_tdof_list_e, phP, B1t, KmP, X1v, B1v);

	HypreBoomerAMG Mpp(KmP);
	Mpp.SetPrintLevel(0);
	
	CGSolver cgPP(MPI_COMM_WORLD);
	cgPP.SetRelTol(1e-6);
	cgPP.SetMaxIter(82);
	cgPP.SetPreconditioner(Mpp);
	cgPP.SetOperator(KmP);
	
	// force Vector
	ParLinearForm Fpt(&fespace);	
	HypreParVector Fpb(&fespace);

	std::unique_ptr<ParLinearForm> Bp2(new ParLinearForm(&fespace));
	Bp2->AddDomainIntegrator(new DomainLFIntegrator(cRp));	
	Bp2->Assemble();		
	Fpt = *Bp2;	// Move the contents of Bp2 into Fpt

	// project values to DBC nodes	
	Kp2->FormLinearSystem(ess_tdof_list_e, phP, Fpt, KmP, X1v, Fpb);

	HypreParVector Xs0(&fespace);


	// ========================================
	// // 	              _     ______ 
	// // 	             | |   |  ____|
	// // 	  _ __   ___ | |_  | |__   
	// // 	 | '_ \ / _ \| __| |  __|  
	// // 	 | |_) | (_) | |_  | |____ 
	// // 	 | .__/ \___/ \__| |______|
	// // 	 | |                       
	// // 	 |_|                       
	// ========================================
	
	double tc1 =(2*t_minus-1.0)/(2*t_minus*(1.0-t_minus));
	double tc2 = 1.0/(2*t_minus*(1.0-t_minus))*Cst1;
	double dffe;

	ParGridFunction phE(&fespace);		// electropot in electrolyte
	ParGridFunction pE0(&fespace);		// values before iteration
		
	ParGridFunction Dmp(&fespace);		// D_minus_plus
	GridFunctionCoefficient cDm(&Dmp);
	
	ParGridFunction kpl(&fespace);		// electrolyte conductivity
	GridFunctionCoefficient cKe(&kpl) ;	
	
	ParGridFunction RpE(&fespace);		// reaction rate for electrolyte
	RpE.Neg();	
	GridFunctionCoefficient cRe(&RpE) ;	

	double BvE = -0.4686;
	phE = BvE;

	// stiffness matrix
	HypreParMatrix Kml;
	
	std::unique_ptr<ParBilinearForm> Kl2(new ParBilinearForm(&fespace)); 		
	Kl2->AddDomainIntegrator(new DiffusionIntegrator(cKe));
	Kl2->Assemble();		

	// assign known values to the DBC nodes	
	ConstantCoefficient dbc_w_Coef(BvE);
	phE.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 		
	Kl2->FormLinearSystem(ess_tdof_list_w, phE, B1t, Kml, X1v, B1v);	

	// // assign known values to the DBC nodes	
	// ConstantCoefficient dbc_e_Coef(BvE);
	// phE.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 		
	// Kl2->FormLinearSystem(ess_tdof_list_e, phE, B1t, Kml, X1v, B1v);	

	HypreBoomerAMG Mpe(Kml);
	Mpe.SetPrintLevel(0);
		
	CGSolver cgPE(MPI_COMM_WORLD);
	cgPE.SetRelTol(1e-6);
	cgPE.SetMaxIter(80);
	cgPE.SetPreconditioner(Mpe);	
	cgPE.SetOperator(Kml);
	
	// force vector
	ParLinearForm Flt(&fespace);
	HypreParVector Flb(&fespace);
	
	
	std::unique_ptr<ParLinearForm> Bl2(new ParLinearForm(&fespace));
	Bl2->AddDomainIntegrator(new DomainLFIntegrator(cRe));	
	Bl2->Assemble();		
	Flt = *Bl2;	// Move the contents of Bl2 into Flt		

	// Kl2->FormLinearSystem(ess_tdof_list_w, phE, Flt, Kml, X1v, Flb);
	Kl2->FormLinearSystem(ess_tdof_list_w, phE, Flt, Kml, X1v, Flb);

							
							
	// Laplace matrix
	HypreParMatrix Kdm;
	
	std::unique_ptr<ParBilinearForm> Kl1(new ParBilinearForm(&fespace));
	Kl1->AddDomainIntegrator(new DiffusionIntegrator(cDm));
	Kl1->Assemble();
	Kl1->FormLinearSystem(boundary_dofs, phE, B1t, Kdm, X1v, B1v);			
	
	ParGridFunction Rxn(&fespace);	

	HypreParVector LpCe(&fespace), Xe0(&fespace);
	HypreParVector RHSl(&fespace);
	
	Rxn = 0.0;
	// Rxn = AvP; 
	// Rxn *= 8.0e-10 * 0.99 *2.25;
	
	double Vcell = BvP - BvE;

	// reaction term
	// ParGridFunction Rxn(&fespace);	
	// Rxn = 8.0e-9;
	 
	// rate constants
	ParGridFunction dPHE(&fespace);		// voltage drop
	ParGridFunction Kfw(&fespace);		// forward reaction constant
	ParGridFunction Kbw(&fespace);		// backward rection constant
	ParGridFunction OCV(&fespace);		// open circuit voltage
	ParGridFunction i0C(&fespace);		// exchange current density	

	OCV = 0.0;
	i0C = 0.0;
	Kfw = 0.0;
	Kbw = 0.0;	

	double errP = 1.0;
	double errE = 1.0;
	int inlp = 0;
	double gErrP, gErrE;
	double dV, sgn;

	int cnt = 0;
	string stri;
	string fname;

	std::cout << "Starting time loop..." << endl;
		
	int indxTick = 0;
	// 	===================================================================
	// 	//   _   _                      _                   _             
	// 	//  | | (_)                    | |                 (_)            
	// 	//  | |_ _ _ __ ___   ___   ___| |_ ___ _ __  _ __  _ _ __   __ _ 
	// 	//  | __| | '_ ` _ \ / _ \ / __| __/ _ \ '_ \| '_ \| | '_ \ / _` |
	// 	//  | |_| | | | | | |  __/ \__ \ ||  __/ |_) | |_) | | | | | (_| |
	// 	//   \__|_|_| |_| |_|\___| |___/\__\___| .__/| .__/|_|_| |_|\__, |
	// 	//                                     | |   | |             __/ |
	// 	//                                     |_|   |_|            |___/ 
	//  ===================================================================   
		
	// timestepping
	int t = 0;
	for (int t = 0; t < 30; t++){
	// 	while ( Vcell > Vcut){
	// while ( Xfr < 0.971 ){
	
	
		// // 	=============================
		// // 	   _____       _____  
		// // 	  / ____|     |  __ \ 		//
		// // 	 | |     _ __ | |__) |
		// // 	 | |    | '_ \|  ___/ 
		// // 	 | |____| | | | |     
		// // 	  \_____|_| |_|_|     
		// // 	==============================
	
		mfem::StopWatch sw;
		sw.Start();
			
		Rxc = Rxn;
		Rxc /= rho;
		// Rxc *= AvP;

		cAp.SetGridFunction(&Rxc);

		Bc2->Update();
		Bc2->Assemble();		
		Fct = *Bc2;	

		// tabulate; mu and mob
		for (int vi = 0; vi < nV; vi++){
			Ctmp = CnP(vi);
			if ( Ctmp < 1.0e-6 ){Ctmp = 1.0e-6;}
			if ( Ctmp > 1.0 ){Ctmp = 1.0;}			
			indxTick = std::floor(Ctmp/d_X);

			if (indxTick < 0) indxTick = 0;
    		if (indxTick > pnt1 - 2) indxTick = pnt1 - 2;

			Mub(vi) = chmPot[indxTick] + (Ctmp-Ticks[indxTick])/d_X * 
				(chmPot[indxTick+1]-chmPot[indxTick]);
			Mob(vi) = Mobil[indxTick] + (Ctmp-Ticks[indxTick])/d_X * 
				(Mobil[indxTick+1]-Mobil[indxTick]);
			Mob(vi) *= psi(vi); // multiply by psi for correct mobility
		}

		// Mob *= psi;		// multiply by psi for correct mobility

		// Mob.Save("Results/mob");
		// Mub.Save("Results/mu");

		// // try
		// GridFunction smoothMu(&fespace);
		// GridFunction smoothMob(&fespace);
		// smoothMu.ProjectGridFunction(Mub);
		// smoothMob.ProjectGridFunction(Mob);

		// cDp.SetGridFunction(&smoothMob);	// update the coefficient for the stiffness matrix

		// Kx2->Update();
		// Kx2->Assemble();
		// Kx2->FormLinearSystem(boundary_dofs, CnP, Fct, KmX, X1v, Fcb);		


		// cDp.SetGridFunction(&Mob);	// update the coefficient for the stiffness matrix

		// vector of CnP	
		CnP.GetTrueDofs(CpV0);
				
		// Lap phi			
		KmX.Mult(CpV0, Lp1);

		Mub.GetTrueDofs(MuV);
		// Mub.GetTrueDofs(MuV); try

		// MuV += Lp1; // mu = mu_b - Lp1	
		MuV += Lp1; // mu = mu_b - Lp1
					
		// K matrix		
		Kc2->Update();
   		Kc2->Assemble();
   		
		
		Kc2->FormLinearSystem(boundary_dofs, Mub, Fct, Kmatp, X1v, Fcb);
		// Kc2->FormLinearSystem(boundary_dofs, Mub, Fct, Kmatp, X1v, Fcb);
   		
		// MuV *= PsVc;  // where PsVc is the nodal vector of psi
		Kmatp.Mult(MuV, Lp2);
		Lp2.Neg();  
		Lp2 *= dt;
				
		// reaction vector
   		Fcb *= dt;

   		Lp2 += Fcb;

		// right hand side
		Mmatp.Mult(CpV0, RHCp);
		RHCp += Lp2;			
	
		// time stepping
		Mp_solver.Mult(RHCp, CpV0) ;

		for (int p = 0; p < nDof; p++) {
			if (PsVc(p) < 1.0e-5) {
				CpV0(p) = Cp0;
			}
			// double C_min = 1.0e-8;
			// double C_max = 0.999999;
			// if (CpV0(p) < 0.0) { CpV0(p) = 0.0; }
			// if (CpV0(p) > 1.0) { CpV0(p) = 1.0; }
		}
		
		// Recover the GridFunction from Vector.
		CnP.Distribute(CpV0);

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

		if (t % 200 == 0 && myid == 0) {
            std::cout << "timestep: " << t
                    << ", Xfr = " << Xfr
                    << std::endl; }

		// cout << "Xrf = " << Xfr << " " << t << endl;

	// }
// 		cout << " -- " << myid << " " << t << " = " << lSum <<  " " << gSum << " " << gtPsi << endl;


		// // ============================
		// //    _____       ______ 
		// //   / ____|     |  ____|
		// //  | |     _ __ | |__   
		// //  | |    | '_ \|  __|  
		// //  | |____| | | | |____ 
		// //   \_____|_| |_|______|
		// // ============================

		Rxe = Rxn;
		Rxe *= (-1.0*t_minus);	

		cAe.SetGridFunction(&Rxe);
		
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

		// std::cout << "L_w = " << L_w << ", infx = " << infx << ", t = " << t << std::endl;

		// need to fix nbcCoef and BCs for Neumann CnE
		// Neumann BC
		ConstantCoefficient nbcCoef(infx);
		ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

		Be2->Assemble();
		Fet = *Be2;
				
		// salt diffusivity in the electrolyte					
		for (int vi = 0; vi < nV; vi++){
			// appendix equation A-21
			De(vi) = pse(vi)*D0*exp(-7.02-830*CnE(vi)+50000*CnE(vi)*CnE(vi));
		}

		cDe.SetGridFunction(&De);	// update the coefficient for the stiffness matrix
		
 		// K matrix		
		Ke2->Update();
   		Ke2->Assemble();
   		Ke2->FormLinearSystem(boundary_dofs, CnE, Fet, Kmate, X1v, Feb);
   		Feb *= dt;		
		
		// Crank-Nicolson matrices	
		TmatR = Add(1.0, Mmate, -0.5*dt, Kmate);		
		TmatL = Add(1.0, Mmate,  0.5*dt, Kmate);		
		
		// vector of CnE				
		CnE.GetTrueDofs(CeV0);		
		
		// solver
		Me_solver.SetOperator(*TmatL);
		
    	TmatR->Mult(CeV0, RHSe);
    	RHSe += Feb;
        	
    	// time stepping
		Me_solver.Mult(RHSe, CeVn) ;
		
		// recover
		CnE.Distribute(CeVn); 
		
		sw.Stop();
		// if (myid == 0) {
		// 	std::cout << "Concentration timestep: " << t 
		// 			<< ", time taken  = " << sw.RealTime() 
		// 			<< " seconds" << std::endl; }


		mfem::StopWatch sw2;
		sw2.Start();

		// check conservation of salt
		if (t%100 == 0){
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
			
			// if (myid==0){cout << CeAvg << " " << L_w << endl;}

		}
		
		sw2.Stop();
		// if (myid == 0) {
		// 	std::cout << "Salt Conservation timestep: " << t
		// 			<< ", time taken  = " << sw2.RealTime() 
		// 			<< " seconds" << std::endl; }
		
		delete TmatR;
		delete TmatL;
		
		
		// ==============================================
		// // 	                      _   _             
		// // 	                     | | (_)            
		// // 	  _ __ ___  __ _  ___| |_ _  ___  _ __  
		// // 	 | '__/ _ \/ _` |/ __| __| |/ _ \| '_ \   //
		// // 	 | | |  __/ (_| | (__| |_| | (_) | | | |
		// // 	 |_|  \___|\__,_|\___|\__|_|\___/|_| |_|
		// ==============================================	
 
		mfem::StopWatch sw3;
		sw3.Start();

		// electrolyte conductivity and RHS	
		for (int vi = 0; vi < nV; vi++){
			dffe = exp(-7.02-830*CnE(vi)+50000*CnE(vi)*CnE(vi));
			Dmp(vi) = pse(vi)*tc1*D0*dffe;
			kpl(vi) = pse(vi)*tc2*D0*dffe*CnE(vi);
		}
		cDm.SetGridFunction(&Dmp);	// update the coefficient for the stiffness matrix

		// Laplace of CnE for the RHS
		Kl1->Update();
		Kl1->Assemble();
		Kl1->FormLinearSystem(boundary_dofs, phE, B1t, Kdm, X1v, B1v);		
	
		// Vector of CnE
		CnE.GetTrueDofs(CeVn) ;
		Kdm.Mult(CeVn, LpCe) ;

		cKe.SetGridFunction(&kpl);	// update the coefficient for the stiffness matrix

		Kl2->Update();
		Kl2->Assemble();	

		// assign known values to the DBC nodes	
		ConstantCoefficient dbc_w_Coef(BvE);
		
		phE.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 		
		Kl2->FormLinearSystem(ess_tdof_list_w, phE, B1t, Kml, X1v, B1v);	
		
		// ConstantCoefficient dbc_e_Coef(BvE);
		// phE.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 		
		// Kl2->FormLinearSystem(ess_tdof_list_e, phE, B1t, Kml, X1v, B1v);		
	
		
		// Solve the system using PCG with hypre's BoomerAMG preconditioner.
		Mpe.SetOperator(Kml);
		cgPE.SetPreconditioner(Mpe);
		cgPE.SetOperator(Kml);

// 		// particle conductivity
// 		// appendix equation A-20
// 		for (int vi = 0; vi < nV; vi++){
// 			kap(vi) = psi(vi)*(0.01929 + 0.7045*tanh(2.399*CnP(vi)) - \
// 				0.7238*tanh(2.412*CnP(vi)) - 4.2106e-6);
// 		}	
// // 		GridFunctionCoefficient cKp(&kap) ;
// 		cKp.SetGridFunction(&kap);
// 		
		// stiffness matrix for phP
// 		std::unique_ptr<ParBilinearForm> Kp2(new ParBilinearForm(&fespace));
// 		Kp2->AddDomainIntegrator(new DiffusionIntegrator(cKp));
		// Kp2->Update();
		// Kp2->Assemble();
// 
		// assign known values to the DBC nodes	
		ConstantCoefficient dbc_e_Coef(BvP);	
				
		// // project values to DBC nodes
		// phP.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 	
		// Kp2->FormLinearSystem(ess_tdof_list_e, phP, B1t, KmP, X1v, B1v);	
		
		// // assign known values to the DBC nodes	
		// // ConstantCoefficient dbc_w_Coef(BvP);	
				
		// // // project values to DBC nodes
		// // phP.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 	
		// // Kp2->FormLinearSystem(ess_tdof_list_w, phP, B1t, KmP, X1v, B1v);

		// // Solve the system using PCG with hypre's BoomerAMG preconditioner.
		// Mpp.SetOperator(KmP);
		// cgPP.SetPreconditioner(Mpp);
		// cgPP.SetOperator(KmP);
		

// 		// rate constants and exchange current density at interface
// 		for (int vi = 0; vi < nV; vi++){
// 			if ( AvB(vi)*dh > 0.0 ){
// 				val = -0.2*(CnP(vi)-0.37)-1.559-0.9376*tanh(8.961*CnP(vi)-3.195);
// 				i0C(vi) = pow(10.0,val)*1.0e-3;
// 				
// 				OCV(vi) = 1.095*CnP(vi)*CnP(vi) - 8.324e-7*exp(14.31*CnP(vi)) + \
// 					4.692*exp(-0.5389*CnP(vi));
// 					
// 				Kfw(vi) = i0C(vi)/(Frd*0.001  )*exp( alp*Cst1*OCV(vi)) ;	
// 				Kbw(vi) = i0C(vi)/(Frd*CnP(vi))*exp(-alp*Cst1*OCV(vi)) ;
// 			}
// 		}	

		// tabulate; OCV and i0
		for (int vi = 0; vi < nV; vi++){
			Ctmp = CnP(vi);
			if ( Ctmp < 1.0e-6 ){Ctmp = 1.0e-6;}
			if ( Ctmp > 9.99999e-1 ){Ctmp = 9.99999e-1;}			
			indxTick = std::floor(Ctmp/d_X);

			if (indxTick < 0) indxTick = 0;
    		if (indxTick > pnt1 - 2) indxTick = pnt1 - 2;

			OCV(vi) = OCVtb[indxTick] + (Ctmp-Ticks[indxTick])/d_X * 
				(OCVtb[indxTick+1]-OCVtb[indxTick]);
			i0C(vi) = i0tb[indxTick] + (Ctmp-Ticks[indxTick])/d_X * 
				(i0tb[indxTick+1]-i0tb[indxTick]);
			i0C(vi) *= 1.0e-3;
			
			Kfw(vi) = i0C(vi)/(Frd*0.001  )*exp( alp*Cst1*OCV(vi)) ;	
			Kbw(vi) = i0C(vi)/(Frd*Ctmp)*exp(-alp*Cst1*OCV(vi)) ;
		}

		// OCV.Save("Results/OCV");
		// i0C.Save("Results/i0");
		// Kfw.Save("Results/Kfw");
		// Kbw.Save("Results/Kbw");

		// convergence residuals
		gErrP = 1.0;
		gErrE = 1.0;
		inlp = 0;

		sw3.Stop();
		// if (myid == 0) {
		// 	std::cout << "Potential timestep: " << t
		// 			<< ", time taken  = " << sw3.RealTime() 
		// 			<< " seconds" << std::endl; }

		// if (t % 100 == 0) {

		mfem::StopWatch sw4;
		sw4.Start();

		// // internal loop
		while (gErrP > 1.0e-8 || gErrE > 1.0e-8 ){

			// Butler-Volmer Equation for Reaction Rate
			for (int vi = 0; vi < nV; vi++){

				if ( AvB(vi)*dh > 0.0 ){
					dPHE(vi) = phP(vi) - phE(vi);
					Rxn(vi) = AvP(vi)*(Kfw(vi)*CnE(vi)*exp(-alp*Cst1*dPHE(vi)) - \
					                   Kbw(vi)*CnP(vi)*exp( alp*Cst1*dPHE(vi)));
				}
				// else {
				// 	Rxn(vi) = 0.0;
				// }
			}

		// 	    std::ofstream outFile("dPHE_values_NOOP.txt"); // open file for writing
		// 		if (!outFile) {
		// 			std::cerr << "Error: Could not open output file.\n";
		// 		} else {
		// 			outFile << std::setprecision(12); // optional for more digits
		// 			outFile << "dPHE values:\n";
		// 			for (int vi = 0; vi < nV; ++vi) {
		// 				outFile << dPHE(vi) << "\n";
		// 			}
		// 			outFile << "\n";
		// 			outFile.close();
		// 		}

			// cout << "min Rxn: " << Rxn.Min() << endl;
			// cout << "max Rxn: " << Rxn.Max() << endl;

			// cout << "min AvB: " << AvB.Min() << endl;
			// cout << "max AvB: " << AvB.Max() << endl;

			// cout << "min dPHE: " << dPHE.Min() << endl;
			// cout << "max dPHE: " << dPHE.Max() << endl;

			// cout << "avg phP: " << phP.Norml2() << endl;
			// cout << "avg phE: " << phE.Norml2() << endl;

			// cout << "min phP: " << phP.Min() << endl;
			// cout << "max phP: " << phP.Max() << endl;
			// cout << "min phE: " << phE.Min() << endl;
			// cout << "max phE: " << phE.Max() << endl;

			// cout << "min CnE: " << CnE.Min() << endl;
			// cout << "max CnE: " << CnE.Max() << endl;

			// cout << "min CnP: " << CnP.Min() << endl;
			// cout << "max CnP: " << CnP.Max() << endl;

			// cout << "min OCV: " << OCV.Min() << endl;
			// cout << "max OCV: " << OCV.Max() << endl;
			// cout << "min i0C: " << i0C.Min() << endl;
			// cout << "max i0C: " << i0C.Max() << endl;
			// cout << "min Kfw: " << Kfw.Min() << endl;
			// cout << "max Kfw: " << Kfw.Max() << endl;
			// cout << "min Kbw: " << Kbw.Min() << endl;
			// cout << "max Kbw: " << Kbw.Max() << endl;


			// ========================================
			// //               _     _____  
			// //              | |   |  __ \ 		//
			// //   _ __   ___ | |_  | |__) |
			// //  | '_ \ / _ \| __| |  ___/ 
			// //  | |_) | (_) | |_  | |     
			// //  | .__/ \___/ \__| |_|     
			// //  | |                       
			// //  |_|                       
			// ========================================		
		
		
			// force vector
			RpP = Rxn;
			RpP *= Frd;	
			cRp.SetGridFunction(&RpP);
		
			Bp2->Assemble();
			Fpt = *Bp2;
			
			// cout << "BvP: " << BvP << endl;


			// project values to DBC nodes
			phP.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 	
			mfem::StopWatch swA, swB;
			swA.Start();
			Kp2->FormLinearSystem(ess_tdof_list_e, phP, Fpt, KmP, X1v, Fpb);	
			swA.Stop();
			if (myid == 0) {
				std::cout << "Kp2 FormLinearSystem time: "
						<< swA.RealTime()
						<< " seconds" << std::endl; }
			
			// // project values to DBC nodes
			// phP.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 	
			// Kp2->FormLinearSystem(ess_tdof_list_w, phP, Fpt, KmP, X1v, Fpb);
			
			pP0 = phP;	
			pP0.GetTrueDofs(Xs0);

			swB.Start();
			cgPP.Mult(Fpb, Xs0);	
			swB.Stop();
			if (myid == 0) {
				std::cout << "cgPP.Mult time: "
						<< swB.RealTime()
						<< " seconds" << std::endl; }
			
			// recover
			phP.Distribute(Xs0);   			

			for (int vi = 0; vi < nV; vi++){
				TmpF(vi) = pow(pP0(vi)-phP(vi),2)*psi(vi);
			}	
				
			errP = 0.0;
			for (int ei = 0; ei < nE; ei++){
				TmpF.GetNodalValues(ei,VtxVal) ;
				val = 0.0;
				for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
				EAvg(ei) = val/nC;						
				errP += EAvg(ei)*EVol(ei) ;
			}
			MPI_Allreduce(&errP, &gErrP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			gErrP /= gtPsi;
			gErrP = pow(gErrP, 0.5);
			
			
	// 		// ========================================
	// 		// // 	              _     ______ 
	// 		// // 	             | |   |  ____|
	// 		// // 	  _ __   ___ | |_  | |__   
	// 		// // 	 | '_ \ / _ \| __| |  __|  
	// 		// // 	 | |_) | (_) | |_  | |____ 
	// 		// // 	 | .__/ \___/ \__| |______|
	// 		// // 	 | |                       
	// 		// // 	 |_|                       
	// 		// ========================================
		
		
			RpE = Rxn;
			RpE.Neg();
			cRe.SetGridFunction(&RpE);	

			// cout << "min RpE: " << RpE.Min() << endl;
			// cout << "max RpE: " << RpE.Max() << endl;

			Bl2->Assemble();		
			Flt = *Bl2;

			phE.ProjectBdrCoefficient(dbc_w_Coef, dbc_w_bdr); 		
			Kl2->FormLinearSystem(ess_tdof_list_w, phE, Flt, Kml, X1v, Flb);

			// cout << "min Flb: " << Flb.Min() << endl;
			// cout << "max Flb: " << Flb.Max() << endl;
		
			// cout << "Essential dofs: " << ess_tdof_list_w.Size() << endl;

			RHSl = Flb;
			RHSl += LpCe;

    		// cout << "norml2 RHSl: " << RHSl.Norml2() << endl;
			// cout << "max RHSl: " << RHSl.Max() << endl;

			pE0 = phE;
			pE0.GetTrueDofs(Xe0);

			cgPE.Mult(RHSl, Xe0);
			
			// cout << "min Xe0: " << Xe0.Min() << endl;
			// cout << "max Xe0: " << Xe0.Max() << endl;
			
			// // recover
			phE.Distribute(Xe0);		

			for (int vi = 0; vi < nV; vi++){
				TmpF(vi) = pow(pE0(vi)-phE(vi),2)*pse(vi);
			}	

			errE = 0.0;	
			for (int ei = 0; ei < nE; ei++){
				TmpF.GetNodalValues(ei,VtxVal) ;
				val = 0.0;
				for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
				EAvg(ei) = val/nC;	
				errE += EAvg(ei)*EVol(ei) ;					
			}	
			MPI_Allreduce(&errE, &gErrE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
			gErrE /= gtPse;
			gErrE = pow(gErrE,0.5);
			
			// inlp += 1;

		} // internal loop

		sw4.Stop();
		if (myid == 0) {
			std::cout << "Advance timestep: " << t
					<< ", time taken  = " << sw4.RealTime() 
					<< " seconds" << std::endl; }

	// // }


		// total reaction current
		sCrnt = 0.0;
		for (int ei = 0; ei < nE; ei++){
			Rxn.GetNodalValues(ei,VtxVal) ;
			val = 0.0;
			for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
			EAvg(ei) = val/nC;
			sCrnt += EAvg(ei)*EVol(ei) ;
		} 			
		MPI_Allreduce(&sCrnt, &gCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



		// dCrnt = abs(gCrnt - gTrgI);
		// if (dCrnt < abs(gTrgI)*0.05){Vsr = 0.025*Vsr0;}
		// else if (dCrnt < abs(gTrgI)*0.10){Vsr = 0.25*Vsr0;}		
		// else {Vsr = 1*Vsr0;}	
		// if (Xfr < 0.051 && gCrnt - gTrgI < 0.0){Vsr = 4*Vsr;}
					
		// adjust BvP for constant C rate loading
		sgn = copysign(1, gTrgI - gCrnt);
		dV = dt*Vsr*sgn;
		// BvP -= dV;
		// phP -= dV;

		BvE += dV;
		phE += dV;
		
		Vcell = BvP - BvE;

// // 		cout << myid << " " << sgn << " " << gTrgI-gCrnt << endl;
		
		tm = tm + dt;

		// std::cout << "timestep: " << t 
		// 		<< ", VCell = " << Vcell 
		// 		<< ", BvP = " << BvP 
		// 		<< ", BvE = " << BvE 
		// 		<< ", gCrnt = " << gCrnt 
		// 		<< std::endl;

	

// // 		if ( t == 5000){
// 		// if ( Xfr > 0.02+cnt*0.005 ){
// 		// 	// output values to a text file
// 		// 	if (myid == 1){
// 		// 		ofstream myfile;
// 		// 		myfile.open ("Output_2D_gra_disk_P01.txt", ios::app);
// 		// 		myfile << cnt << "  " << t << "  " << tm << "  " << Xfr << "  " << Vcell << "  " \
// 		// 			<< "  " << gCrnt << "  " << gTrgI << "  " << CeAvg << endl;
// 		// 		myfile.close();
// 		// 	}
			
// 		// 	stri = to_string(1000+cnt);
		
// 		// 	fname = "output1/CnP_96x96_P01_" + stri;			
// 		// 	const char *test1 = fname.c_str();
// 		// 	CnP.Save(test1);

// 		// 	fname = "output1/CnE_96x96_P01_" + stri;			
// 		// 	const char *test2 = fname.c_str();
// 		// 	CnE.Save(test2);
			
// 		// 	fname = "output1/phP_96x96_P01_" + stri;			
// 		// 	const char *test3 = fname.c_str();
// 		// 	phP.Save(test3);			

// 		// 	fname = "output1/phE_96x96_P01_" + stri;			
// 		// 	const char *test4 = fname.c_str();
// 		// 	phE.Save(test4);
						
// 		// 	cnt += 1;
			
// 		// 	MPI_Barrier(MPI_COMM_WORLD);
// 		// }
		
// 		// t += 1;
		

// // 		// if ( Vcell < Vcut ){
// // 		// 	// output values to a text file
// // 		// 	if (myid == 1){
// // 		// 		ofstream myfile;
// // 		// 		myfile.open ("Output_2D_gra_disk_P01.txt", ios::app);
// // 		// 		myfile << cnt << "  " << t << "  " << tm << "  " << Xfr << "  " << Vcell << "  " \
// // 		// 			<< "  " << gCrnt << "  " << gTrgI << "  " << CeAvg << endl;
// // 		// 		myfile.close();
// // 		// 	}
			
// // 		// 	stri = to_string(1000+cnt);
 		
// // 		// 	fname = "output1/CnP_96x96_P02_" + stri;			
// // 		// 	const char *test1 = fname.c_str();
// // 		// 	CnP.Save(test1);

// // 		// 	fname = "Output1/CnE_80x96_S02_" + stri;			
// // 		// 	const char *test2 = fname.c_str();
// // 		// 	CnE.Save(test2);
			
// // 		// 	fname = "Output1/phP_80x96_S02_" + stri;			
// // 		// 	const char *test3 = fname.c_str();
// // 		// 	phP.Save(test3);			

// // 		// 	fname = "Output1/phE_80x96_S02_" + stri;			
// // 		// 	const char *test4 = fname.c_str();
// // 		// 	phE.Save(test4);
						
// // 		// 	cnt += 1;				
// // 		// }
		
		
// // 	// if (t%50 == 1 && myid == 1){cout << t << " - " << Xfr << " - " << tm << " - " << \
// // 	// 	Vcell << " xx " << inlp-1 << " - " << gErrP << " - " << gErrE << \
// // 	// 	" " << gCrnt << " -- " << gTrgI << endl;}
// // // 	if (myid == 1 ){cout << t << "  " << Xfr << "  " << tm << "  " << Vcell << endl;}
		
		if (t % 1 == 0 && myid == 0) {
            std::cout << "timestep: " << t
                    << ", Xfr = " << Xfr
                    << ", VCell = " << Vcell << ", BvE = " << BvE
                    << std::endl;
        }



	} // time iteration loop


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

    CnP.Save("Results/CnCH");
	
    CnE.Save("Results/CnE");
    phP.Save("Results/phP");
    phE.Save("Results/phE");    

    phP *= psi;
    CnE *= pse;
    phE *= pse;
	CnP *= psi;

    CnP.Save("Results/pCnCH");
    CnE.Save("Results/pCnE");
    phP.Save("Results/pphP");
    phE.Save("Results/pphE");
    Rxn.Save("Results/Rxn");

	ParGridFunction muField(&fespace);
	muField.SetFromTrueDofs(MuV);  // transfer MuV to nodal field
	muField.Save("Results/muCH");  
	Mob.Save("Results/MobCH");

}
	// Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

	// End timing and output the total program execution time
    auto program_end = high_resolution_clock::now();
    std::cout << "Total Program Time: " 
              << duration_cast<seconds>(program_end - program_start).count() 
              << " seconds" << std::endl;
	
	
	
   return 0;
// }	

}