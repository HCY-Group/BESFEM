#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include "TimeDepOpers.hpp"
#include "readtiff.h"

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
	int nz = tiffdata.size();
	int ny = tiffdata[0].size();
	int nx = tiffdata[0][0].size();
	double sx = nx;  //make dx = 1
	double sy = ny;  //make dy = 1
	bool generate_edges = false;
	bool sfc_ordering = false;
	Mesh gmesh = Mesh::MakeCartesian2D(nx-1, ny-1, Element::QUADRILATERAL, generate_edges, sx, sy, sfc_ordering);
	gmesh.EnsureNCMesh(true);

	// Create global FE space for Voxel Data
	int order = 1;
	H1_FECollection gFec(order, gmesh.Dimension());
	FiniteElementSpace gFespace(&gmesh, &gFec);

	// global grid function for voxel data
	cout << "Defining Voxel GridFunction" << endl;
	GridFunction gVox(&gFespace);
	for(int j=0; j<ny; j++){
		for(int i=0; i<nx; i++){
			int idx = i+nx*j;
			//cout << "idx = " << i+ny*j << endl;
			//cout << "gVox[idx] = " << gVox[idx] << endl;
			//cout << "tiffdata[i][j][0] = " << tiffdata[i][j][0] << endl;
			gVox[idx] = tiffdata[0][j][i];
		}
	}
	
	// double check our read-in: output gVox to Paraview
	cout << "PRINTING gVox" << endl;
	ParaViewDataCollection *pd = NULL;
	pd = new ParaViewDataCollection("gVoxelData", &gmesh);
	pd->RegisterField("gVox", &gVox);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
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
	for (int t = 0; t < 30; t++){
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

	// Concentration GridFunction
	ParGridFunction Cn(&fespace);
	Cn = 0.0;
	// Dirichlet boundary value
	double Bv = 1.0;
	
	// Define order parameter psi from distance function
	// This will use diffusion, so here we use H1 Finite Element Space again
	ParGridFunction psi(&fespace);
	psi.ProjectGridFunction(d);
	double zeta = 0.75;
	for (int vi = 0; vi < nV; vi++){
		psi(vi) = 0.5*( 1.0+tanh(psi(vi)/zeta) );
	}
	
	// Indicate Dirichlet boundary conditions
	Array<int> dbc_bdr(pmesh.bdr_attributes.Max());
	dbc_bdr = 0; dbc_bdr[2] = 1;
	// Node labels of Dirichlet BC
	Array<int> ess_tdof_list(0);
	fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);
	// Assign Dirichlet values
	ConstantCoefficient dbc_Coef(Bv);
	Cn.ProjectBdrCoefficient(dbc_Coef, dbc_bdr);
	
	// stiffness matrix and forcing vector
	std::unique_ptr<ParBilinearForm> KCn(new ParBilinearForm(&fespace));
	KCn->AddDomainIntegrator(new DiffusionIntegrator);
	KCn->Assemble();
	
	ParLinearForm FctCn(&fespace);
	HypreParMatrix KmatCn;
	HypreParVector X1vCn(&fespace);
	HypreParVector FcbCn(&fespace);
	KCn->FormLinearSystem(ess_tdof_list, Cn, FctCn, KmatCn, X1vCn, FcbCn);
	cout << "b Size: " << FctCn.Size() << endl;
	cout << "B Size: " << FcbCn.Size() << endl;

	ParGridFunction psi0(&fespace);
	KCn->FormLinearSystem(ess_tdof_list, Cn, psi, KmatCn, X1vCn, psi0);
	
	// Initialize TimeDependentOperator and ODESolver
	ConductionOperator operCon(psi0, KmatCn, FcbCn, ess_tdof_list);
	//ConductionOperator operCon(psi, KmatCn, FcbCn);
	
	//ODESolver *ode_solver_Con = new ForwardEulerSolver;
	ODESolver *ode_solver_Con = new BackwardEulerSolver;
	ode_solver_Con->Init(operCon);
	
	// Time Step
	HypreParVector Cn0(&fespace);
	t_ode = 0.0;
	dt = 0.0005;
	for (int t = 0; t < 100; t++){
		Cn.GetTrueDofs(Cn0);
		ode_solver_Con->Step(Cn0, t_ode, dt);
		// update only in the region of interest
		for (int vi = 0; vi < nV; vi++){
			if (psi0(vi) < 1.0e-5){
				Cn0(vi) = 0.0;
			}
		}
		Cn.Distribute(Cn0);
		cout << "iter: " << t << " max Cn: " << Cn0.Max() << endl;
		cout << "iter: " << t << " min Cn: " << Cn0.Min() << endl;
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

	return 0;

}
