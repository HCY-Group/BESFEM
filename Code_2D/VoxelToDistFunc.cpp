#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include "TimeDepOpers.hpp"
#include "readtiff.h"
#include "MeshMaker.hpp"
#include "VoxelSolver.hpp"

using namespace std;
using namespace mfem;





int main(int argc, char *argv[])
{
	// Initialize MPI and HYPRE
	Mpi::Init(argc, argv);
	Hypre::Init();
	
	// ======================================
	// READ IN TIFF FILE
	// ======================================
	cout << "LOADING IN TIFF" << endl;
	const char* tiffname ="II_1_bin.tif";
	Constraints args;
	args.Depth_begin = 0;	//only read in one slice for 2D data
	args.Depth_end = 1;	//only read in one slice for 2D data
	// get a smaller subset so it runs faster
	args.Row_begin    = 0;
	args.Row_end      = 80;
	args.Column_begin = 0;
	args.Column_end   = 120;
	TIFFReader reader(tiffname,args);
	reader.readinfo();
	std::vector<std::vector<std::vector<int>>> tiffdata;
	tiffdata = reader.getImageData();
	
	cout << "tiff size: "            << tiffdata.size()         << endl;	
	cout << "tiff[0] size: "         << tiffdata[0].size()      << endl;
	cout << "tiff[0][0] size: "      << tiffdata[0][0].size()   << endl;
	cout << "tiffdata[0][0][0] = "   << tiffdata[0][0][0]       << endl;
	/*for (int i = 0; i < tiffdata.size(); i++) {
		for (int j = 0; j < tiffdata[i].size(); j++) {
			for (int k = 0; k < tiffdata[i][j].size(); k++) {
				cout << tiffdata[i][j][k] << endl;
			}
		}
	}
	*/


	// ======================================
	// SERIAL MESH AND GLOBAL DATA
	// ======================================
	cout << "MAKING MESH" << endl;
	MeshMaker maker(tiffdata);
	maker.MakeGlobalMesh();
	maker.Make_H1_FESpace();
	//maker.AssignGlobalValues();
	
	VoxelSolver solver(maker.GetGlobalFESpace());
	solver.AssignGlobalValues(tiffdata);
	int order = 1;

	// weirdness to get the paraview save to work.
	// Can't just return the grid function from the class
	// we have to create the FiniteElementSpace again from scratch, starting with the FECollection
	cout << "returning gVox" << endl;
	
	//GridFunction gVox3(*maker.GetGlobalVox());
	GridFunction gVox3(*solver.GetGlobalVox());
	
	//GridFunction gVox2(maker.GetGlobalVox()->FESpace());
	//cout << "HERE" << maker.GetGlobalVox()->FESpace()->GetMesh() << endl;;
	
	
	H1_FECollection gFec2(order, maker.GetGlobalMesh()->Dimension());
	FiniteElementSpace gFespace2(maker.GetGlobalMesh(), &gFec2);
	GridFunction gVox2(&gFespace2);
	
	//FiniteElementSpace *gFespace3 = maker.GetGlobalVox()->FESpace();
	//cout << "HERE" << gFespace3 << endl;
	//FiniteElementSpace gFespace2(*gFespace3);
	//FiniteElementSpace gFespace2(*maker.GetGlobalVox()->FESpace());
	//GridFunction gVox2(&gFespace);
	//GridFunction gVox2(&gFespace2);
	//gVox2 = *maker.GetGlobalVox();

	//for (int idx=0; idx<gVox2.Size(); idx++){
	//	gVox2[idx] = gVox3[idx];
	//	cout << idx << " " << gVox2[idx] << endl;
	//}
	gVox2 = gVox3;
	
	cout << "returning Mesh" << endl;
	Mesh gmesh2(*maker.GetGlobalMesh());

	// double check our read-in: output gVox to Paraview
	cout << "PRINTING gVox" << endl;
	cout << gVox2.Size() << endl;
	cout << gVox2.OwnsData() << endl;
	cout << gVox2[1] << endl;
	cout << gmesh2.Dimension() << endl;
	cout << gmesh2.GetNV() << endl;

	ParaViewDataCollection *pd = NULL;
	cout << "before" << endl;
	//pd = new ParaViewDataCollection("gVoxelData", &gmesh2);
	pd = new ParaViewDataCollection("gVoxelData", maker.GetGlobalMesh());
	//pd->RegisterField("gVox", &gVox);
	pd->RegisterField("gVox", &gVox2);
	//pd->RegisterField("gVox", maker.GetGlobalVox());
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	cout << "after" << endl;
	pd->Save();
	cout << "after" << endl;
	delete pd;
		
	/*
	// glvis -m gmesh.mesh -g gvox.gf
	gVox.Save("gvox.gf");
	gmesh.Save("gmesh.mesh");
	*/

	/*
	char vishost[] = "localhost";
 	int  visport   = 19916;
      	socketstream sol_sock(vishost, visport);
      	sol_sock.precision(8);
      	sol_sock << "solution\n" << gmesh << gVox << flush;
	*/
	/*

	char vishost[] = "localhost";
 	int  visport   = 19916;
      	socketstream sol_sock(vishost, visport);
      	sol_sock << "parallel " << Mpi::WorldSize() << " " << Mpi::WorldRank() << "\n";
      	sol_sock.precision(8);
      	sol_sock << "solution\n" << gmesh << gVox << flush;
	*/











	// ======================================
	// LOCAL (PARALLEL) MESH AND DATA
	// ======================================
	// make variable names match from up above...
	Mesh gmesh(*maker.GetGlobalMesh());
	GridFunction gVox(gVox2);
	
	cout << "starting parallel mesh" << endl;
	ParMesh pmesh(MPI_COMM_WORLD, gmesh);
	int nV = pmesh.GetNV();			//number of vertices
	int nE = pmesh.GetNE();			//number of elements
	int nC = pow(2,pmesh.Dimension());	//number of corner vertices
	
	// Define finite element space
	H1_FECollection fec(order, pmesh.Dimension());
	ParFiniteElementSpace fespace(&pmesh, &fec);

	// Map local to global element indices
	Array<HYPRE_BigInt> E_L2G;
	pmesh.GetGlobalElementIndices(E_L2G);

	// Local GridFunction
	cout << "DEFINING LOCAL GRIDFUNCTION" << endl;
	ParGridFunction Vox(&fespace);
	
	Array<int> gVTX(nC);	//global indices of corner vertices
	Array<int> VTX(nC);	//local indices of corner vertices
	int gei;			//global element indices
	int ei;			//local element indices
	for (ei=0; ei<nE; ei++){
		gei = E_L2G[ei];

		gmesh.GetElementVertices(gei,gVTX);
		pmesh.GetElementVertices(ei,VTX);
	
		for (int vi = 0; vi<nC; vi++){
			Vox(VTX[vi]) = gVox(gVTX[vi]);
		}
	}
	


	// double check our read-in: output Vox to Paraview
	cout << "PRINTING OUT Vox" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("pVoxelData", &pmesh);
	pd->RegisterField("Vox", &Vox);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;






	// ======================================
	// SMOOTH THE DATA USING ALLEN-CAHN
	// ======================================
	cout << "SMOOTHING WITH ALLEN-CAHN" << endl;
	
	// Natural (Neumann) boundary conditions
	Array<int> boundary_dofs;
	
	ParLinearForm Fct(&fespace);
	HypreParVector Fcb(&fespace);
	HypreParVector X1v(&fespace);
	
	// mass matrix
	ParGridFunction ones(&fespace);	
	ones = 1.0;

	// stiffness matrix
	HypreParMatrix Kmat;
	ParGridFunction Mob(&fespace);
	Mob = 0.64;
	GridFunctionCoefficient cMob(&Mob);
	std::unique_ptr<ParBilinearForm> K(new ParBilinearForm(&fespace));
	K->AddDomainIntegrator(new DiffusionIntegrator(cMob));
	K->Assemble();
	K->FormLinearSystem(boundary_dofs, Vox, Fct, Kmat, X1v, Fcb);
	cout << "b max: " << Fct.Max() << endl;
	cout << "B max: " << Fcb.Max() << endl;

	// TimeDependentOperator and ODESolver
	ConductionOperator oper(ones, Kmat, Fcb);
	ODESolver *ode_solver = new ForwardEulerSolver;
	//ODESolver *ode_solver = new BackwardEulerSolver;
	ode_solver->Init(oper);
	
	// time step
	ParGridFunction Pot(&fespace);
	HypreParVector Vox0(&fespace);
	double t_ode = 0.0;
	double dt = 0.05;
	for (int t = 0; t < 20; t++){
		
		//forcing function
		for (int vi = 0; vi < nV; vi++){
			Pot(vi) = 2.0*Vox(vi)*(1.0-Vox(vi))*(1.0-2.0*Vox(vi));
		}
		GridFunctionCoefficient cPot(&Pot);
		std::unique_ptr<ParLinearForm> Bc(new ParLinearForm(&fespace));
		Bc->AddDomainIntegrator(new DomainLFIntegrator(cPot));
		Bc->Assemble();
		Fct = std::move(*Bc);

		//Update Parameters and Solve
		Vox.GetTrueDofs(Vox0);
		K->FormLinearSystem(boundary_dofs, Vox, Fct, Kmat, X1v, Fcb);
		oper.UpdateParams(Kmat, Fcb);
		ode_solver->Step(Vox0, t_ode, dt);
		Vox.Distribute(Vox0);
		
	}
	//oper.~ConductionOperator();
	
	// Output Vox to Paraview
	cout << "PRINTING OUT smoothed Vox" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("SmoothVox", &pmesh);
	pd->RegisterField("Vox", &Vox);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;



	// ======================================
	// FIND DISTANCE FUNCTION BY REINITIALIZING LEVEL SET
	// ======================================
	cout << "FINDING DISTANCE FUNCTION WITH LEVEL SET REINITIALIZATION" << endl;
	
	// Need to use Discontinuous Galerkin (DG) for upwinding (look at MFEM examples 9 and 18)
	DG_FECollection fec_dg(order, pmesh.Dimension(), BasisType::GaussLobatto);
	ParFiniteElementSpace fespace_dg(&pmesh, &fec_dg);
	ParFiniteElementSpace dfespace_dg(&pmesh, &fec_dg, pmesh.Dimension(), Ordering::byNODES); //X1X2X3.....,Y1Y2Y3.....,Z1Z2Z3......

	// Define new grid functions
	ParGridFunction d(&fespace_dg);
	ParGridFunction c(&dfespace_dg);
	ParGridFunction cx(&fespace_dg);
	ParGridFunction cy(&fespace_dg);
	cout << "size of d: " << d.Size() << endl;
	cout << "size of c: " << c.Size() << endl;
	cout << "size of Vox: " << Vox.Size() << endl;
	
	// Set d equal to Vox for initial conditions.
	// Since they are on different FESpaces, we use ProjectGridFunction()
	d.ProjectGridFunction(Vox);
	// Sign function
	ParGridFunction sgn(&fespace_dg);
	for (int vi = 0; vi < sgn.Size(); vi++){
		sgn(vi) = (d(vi)<0.5) ? -1 : 1;
	}
	
	HypreParVector *D = d.GetTrueDofs();
	// Time Stepping
	for (int t = 0; t < 10; t++){
		//gradient of d
		ParGridFunction gdX(&fespace_dg);
		ParGridFunction gdY(&fespace_dg);
		d.GetDerivative(1,0,gdX);
		d.GetDerivative(1,1,gdY);
		
		//calculate c ("velocity")
		ParGridFunction mgGd(&fespace_dg);
		for (int vi = 0; vi < mgGd.Size(); vi++){
			// magnitude of gradient
			mgGd(vi) = sqrt( gdX(vi)*gdX(vi) + gdY(vi)*gdY(vi) );
			// c_x
			c(vi) = sgn(vi)*gdX(vi)/mgGd(vi);
			cx(vi) = c(vi);
			// c_y
			c(vi+mgGd.Size()) = sgn(vi)*gdY(vi)/mgGd(vi);
			cy(vi) = c(vi+mgGd.Size());
		}
		
		//calculate M and K matrices and b vector
		real_t alpha = -1.0;
		ParBilinearForm *m = new ParBilinearForm(&fespace_dg);
		m->AddDomainIntegrator(new MassIntegrator);
		
		ParBilinearForm *k = new ParBilinearForm(&fespace_dg);
		VectorGridFunctionCoefficient cCoef(&c);
   		k->AddDomainIntegrator(new ConvectionIntegrator(cCoef, alpha));
  		k->AddInteriorFaceIntegrator(
      			new NonconservativeDGTraceIntegrator(cCoef, alpha));
   		//k->AddBdrFaceIntegrator(
      		//	new NonconservativeDGTraceIntegrator(cCoef, alpha));
		
		//Add a (small) diffusive term to allow Neumann BCs to be enforced weakly by FEM
		ConstantCoefficient diffCoef(1e-8);
		k->AddDomainIntegrator(new DiffusionIntegrator(diffCoef));
		k->AddBdrFaceIntegrator(
			new DGDiffusionIntegrator(diffCoef,-1,-1));
		k->AddInteriorFaceIntegrator(
			new DGDiffusionIntegrator(diffCoef,-1,-1));
		//k->AddBdrFaceIntegrator(
		//	new DiffusionIntegrator);
		//k->AddBoundaryIntegrator(
		//	new DiffusionIntegrator);

		/*
		cx.Neg();
		GridFunctionCoefficient c_cx(&cx);
		cy.Neg();
		GridFunctionCoefficient c_cy(&cy);
		k->AddBdrFaceIntegrator(
			new BoundaryMassIntegrator(c_cx));
		*/

		ParLinearForm *b = new ParLinearForm(&fespace_dg);
		GridFunctionCoefficient sgnCoef(&sgn);
		b->AddDomainIntegrator(new DomainLFIntegrator(sgnCoef));
		
		//ParGridFunction zeros(&fespace_dg);
		//zeros = 5.0;
		//GridFunctionCoefficient inflow(&cy);
		//b->AddBdrFaceIntegrator(
		//	new BoundaryFlowIntegrator(inflow, cCoef, alpha));
		
		/*
		Array<int> ewbdr(pmesh.bdr_attributes.Max());
		ewbdr = 0; ewbdr[0]=1; ewbdr[2]=1;
		Array<int> nsbdr(pmesh.bdr_attributes.Max());
		nsbdr = 0; nsbdr[1]=1; nsbdr[3]=1;
		//cx.Neg();
		GridFunctionCoefficient c_cx(&cx);
		//cy.Neg();
		GridFunctionCoefficient c_cy(&cy);
		b->AddBdrFaceIntegrator(
			new BoundaryFlowIntegrator(c_cx, cCoef, alpha), ewbdr);
		b->AddBdrFaceIntegrator(
			new BoundaryFlowIntegrator(c_cy, cCoef, alpha), nsbdr);
		//b->AddBdrFaceIntegrator(new BoundaryLFIntegrator(inflow));
		*/

   		int skip_zeros = 0;
   		m->Assemble();
   		k->Assemble(skip_zeros);
   		b->Assemble();
   		m->Finalize();
   		k->Finalize(skip_zeros);

   		HypreParVector *B = b->ParallelAssemble();
		
		// TimeDependentOperator and ODESolver
		FE_Evolution adv(*m, *k, *B);
		double t_ode = 0.0;
		double dt = 0.1;
		ODESolver *ode_solver_dg = new ForwardEulerSolver;
		ode_solver_dg->Init(adv);
		ode_solver_dg->Step(*D, t_ode, dt);
		cout << "iter: " << t << " max distance: " << D->Max() << endl;
		cout << "iter: " << t << " min distance: " << D->Min() << endl;
		
		//free memory
		delete m;
		delete k;
		delete b;
		delete B;
		delete ode_solver_dg;
	
	}
	d.Distribute(D);
	
	
	// Output Distance to Paraview
	cout << "PRINTING OUT DistanceFunction" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("DstFun", &pmesh);
	pd->RegisterField("Dst", &d);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;

	// Output Advection Velocity to Paraview
	//cout << "PRINTING OUT DistanceFunction" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("AdvVel", &pmesh);
	pd->RegisterField("Vel", &c);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;





	// ======================================
	// FIND CONNECTIVITY
	// ======================================

	// Use d to define psi
	ParGridFunction psi(&fespace);
	psi.ProjectGridFunction(d);
	psi -= 0.5; // Center about 0
	for (int vi = 0; vi < nV; vi++){
		psi(vi) = 0.5*( tanh(psi(vi)) + 1.0 );
	}

	// Output psi to Paraview
	cout << "PRINTING OUT Electrolyte Concentration" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("psi", &pmesh);
	pd->RegisterField("psi", &psi);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;

	for (int ConIter = 0; ConIter < 2; ConIter++){
		if (ConIter==1){
			psi.Neg();
			psi += 1.0;
		}
		// Concentration
		ParGridFunction Cn(&fespace);
		Cn = 0.0;

		// Dirichlet Boundary conditions
		Array<int> dbc_bdr(pmesh.bdr_attributes.Max());
		dbc_bdr = 0; 
		if (pmesh.Dimension()==2) {
			if (ConIter==0){
				dbc_bdr[2] = 1; //north?
			} else {
				dbc_bdr[0] = 1; //south?
			}
		} else {
			if (ConIter==0){
				dbc_bdr[3] = 1; //north?
			} else {
				dbc_bdr[1] = 1; //south?
			}
		}	
		Array<int> ess_tdof_list(0);
		fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);
		//ConstantCoefficient dbcCoef(1.0);
		ParGridFunction dbcval(&fespace);
		dbcval = 0.0;
		for (int vi = 0; vi < nV; vi++){
			if ( psi(vi) > 0.6 ){
				dbcval(vi) = 1.0;
			}
		}
		GridFunctionCoefficient dbcCoef(&dbcval);
		Cn.ProjectBdrCoefficient(dbcCoef, dbc_bdr);
		
		// Stiffness matrix
		GridFunctionCoefficient psiCoef(&psi);
		ParBilinearForm KCn(&fespace);
		KCn.AddDomainIntegrator(new DiffusionIntegrator(psiCoef));
		KCn.Assemble();

		// Linear form b
		ParLinearForm b(&fespace);
		b.Assemble();

		// Form Linear System
		HypreParMatrix KmatCn;
		HypreParVector X(&fespace);
		HypreParVector B(&fespace);
		KCn.FormLinearSystem(ess_tdof_list, Cn, b, KmatCn, X, B);

		// TimeDependentOperator and ODESolver
		ConductionOperator operCn(psi, KmatCn, B, ess_tdof_list);
		ODESolver *ode_solverCn = new BackwardEulerSolver;
		ode_solverCn->Init(operCn);

		t_ode = 0.0;
		dt = 0.05;
		for (int t = 0; t < 100; t++){
			ode_solverCn->Step(X, t_ode, dt);
			
			// Accelerate the diffusion
			KCn.RecoverFEMSolution(X, b, Cn);  //Cast X back to Cn
			for (int vi = 0; vi < nV; vi++){
				if ( Cn(vi) > 1.0e-2 && psi(vi) > 0.6 ){
					Cn(vi) = 1.0;  //Modify Cn
				}
			}
			Cn.ProjectBdrCoefficient(dbcCoef, dbc_bdr); //Reapply Dirichlet BCs
			KCn.FormLinearSystem(ess_tdof_list, Cn, b, KmatCn, X, B, 1); //Cast Cn back onto X, making sure to copy interior
			
				
			cout << "Max Cn: " << X.Max() << endl;
			cout << "Min Cn: " << X.Min() << endl;
		}
		KCn.RecoverFEMSolution(X, b, Cn);

		// Output Cn to Paraview
		cout << "PRINTING OUT Electrolyte Concentration" << endl;
		//ParaViewDataCollection *pd = NULL;
		pd = NULL;
		if (ConIter==0){
			pd = new ParaViewDataCollection("Conc_p", &pmesh);
			pd->RegisterField("Cn_p", &Cn);
		} else {
			pd = new ParaViewDataCollection("Conc_e", &pmesh);
			pd->RegisterField("Cn_e", &Cn);
		}
		pd->RegisterField("psi", &psi);
		pd->SetLevelsOfDetail(order);
		pd->SetDataFormat(VTKFormat::BINARY);
		pd->SetHighOrderOutput(true);
		pd->SetCycle(0);
		pd->SetTime(0.0);
		pd->Save();
		delete pd;
	}


















	/*
	// Concentration GridFunction
	ParGridFunction Cn(&fespace);
	Cn = 0.0;
	
	// Define order parameter psi from distance function
	// This will use diffusion, so here we use H1 Finite Element Space again
	ParGridFunction psi(&fespace);
	psi.ProjectGridFunction(d);
	double zeta = 0.75;
	for (int vi = 0; vi < nV; vi++){
		psi(vi) = 0.5*( 1.0+tanh(psi(vi)/zeta) );
	}
	psi = 1.0;
	cout << "Min psi: " << psi.Min() << endl;
	cout << "Max psi: " << psi.Max() << endl;
	
	// Indicate Dirichlet boundary conditions
	Array<int> dbc_bdr(pmesh.bdr_attributes.Max());
	dbc_bdr = 0; dbc_bdr[2] = 1;
	// Node labels of Dirichlet BC
	Array<int> ess_tdof_list(0);
	fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);
	// Assign Dirichlet values
	ConstantCoefficient dbc_Coef(1.0);
	Cn.ProjectBdrCoefficient(dbc_Coef, dbc_bdr);
	
	// stiffness matrix and forcing vector
	GridFunctionCoefficient psi_coef(&psi);
	std::unique_ptr<ParBilinearForm> KCn(new ParBilinearForm(&fespace));
	KCn->AddDomainIntegrator(new DiffusionIntegrator(psi_coef));
	KCn->Assemble();
	
	ParLinearForm FctCn(&fespace);
	HypreParMatrix KmatCn;
	HypreParVector X1vCn(&fespace);
	HypreParVector FcbCn(&fespace);
	//KCn->FormLinearSystem(ess_tdof_list, Cn, FctCn, KmatCn, X1vCn, FcbCn);
	KCn->FormLinearSystem(boundary_dofs, Cn, FctCn, KmatCn, X1vCn, FcbCn);
	cout << "b max: " << FctCn.Max() << endl;
	cout << "B max: " << FcbCn.Max() << endl;

	// Initialize TimeDependentOperator and ODESolver
	ConductionOperator operCon(psi, KmatCn, FcbCn, ess_tdof_list);
	//ConductionOperator operCon(psi, KmatCn, FcbCn, boundary_dofs);
	//operCon.~ConductionOperator();
	
	ODESolver *ode_solver_Con = new ForwardEulerSolver;
	//ODESolver *ode_solver_Con = new BackwardEulerSolver;
	ode_solver_Con->Init(operCon);
	
	// Output Cn to Paraview
	cout << "PRINTING OUT Electrolyte Concentration" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("InitConn_E", &pmesh);
	pd->RegisterField("Conn_E", &Cn);
	pd->RegisterField("psi", &psi);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;

	// Time Step
	HypreParVector Cn0(&fespace);
	t_ode = 0.0;
	dt = 0.0005;
	for (int t = 0; t < 1; t++){
		//Cn.GetTrueDofs(Cn0);
		//ode_solver_Con->Step(Cn0, t_ode, dt);
		//Cn.Distribute(Cn0);
		
		ode_solver_Con->Step(X1vCn, t_ode, dt);
		KCn->RecoverFEMSolution(X1vCn, FctCn, Cn);
		
		// update only in the region of interest
		//for (int vi = 0; vi < nV; vi++){
		//	if (psi(vi) < 1.0e-5){
		//		Cn(vi) = 0.0;
		//	}
		//}
		
		cout << "iter: " << t << " max Cn: " << Cn.Max() << endl;
		cout << "iter: " << t << " min Cn: " << Cn.Min() << endl;
	}
	
	// Output Cn to Paraview
	cout << "PRINTING OUT Electrolyte Concentration" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("Conn_E", &pmesh);
	pd->RegisterField("Conn_E", &Cn);
	pd->RegisterField("psi", &psi);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;

	cout << "Here at end" << endl;
	*/

	return 0;

}
