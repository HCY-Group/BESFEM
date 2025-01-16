#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include "TimeDepOpers.hpp"
#include "readtiff.h"
#include "MeshMaker.hpp"
#include "VoxelSolver.hpp"
#include "VoxelSolver_DG.hpp"

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
	maker.MakeParallelMesh();
	maker.Make_H1_FESpace_Parallel();
	maker.TestGetFE();
	
	
	//VoxelSolver solver(maker.GetGlobalFESpace());
	VoxelSolver solver(maker.GetGlobalFESpace(), maker.GetParallelFESpace());
	solver.AssignGlobalValues(tiffdata);
	solver.ParaviewSave("gVoxelData","gVox",solver.GetGlobalVox());
	
	






	// ======================================
	// LOCAL (PARALLEL) MESH AND DATA
	// ======================================
	solver.MapGlobalToLocal(maker.GetGlobalMesh(),maker.GetParallelMesh());
	solver.ParaviewSave("pVoxelData","Vox",solver.GetParallelVox());
	




	// ======================================
	// SMOOTH THE DATA USING ALLEN-CAHN
	// ======================================
	cout << "SMOOTHING WITH ALLEN-CAHN" << endl;
	
	// Natural (Neumann) boundary conditions
	Array<int> boundary_dofs;
	
	// Domain Parameter
	//ParGridFunction ones(&fespace);	
	ParGridFunction ones(maker.GetParallelFESpace());
	ones = 1.0;
	
	// Diffusion Coefficient
	//ParGridFunction Mob(&fespace);
	ParGridFunction Mob(maker.GetParallelFESpace());
	Mob = 0.64;
	

	solver.InitMatricesAndTimeDepOpers(boundary_dofs, Mob, ones);
	
	// time step
	double t_ode = 0.0;
	double dt = 0.05;
	for (int t = 0; t < 20; t++){

		solver.UpdateLinearForm_DoubleWellPotential();
		
		solver.UpdateSystemAndSolve(boundary_dofs, t_ode, dt);
	}
	//oper.~ConductionOperator();
	
	// Output Vox to Paraview
	solver.ParaviewSave("SmoothVox","Vox",solver.GetGlobalVox());













	// make variable names match from up above...
	ParMesh pmesh(*maker.GetParallelMesh());
	//ParGridFunction Vox2(*solver.GetParallelVox());
	ParGridFunction Vox(*solver.GetParallelVox());
	// TODO: When we construct fespace, do we need to also copy ParMesh and FiniteElementCollection???
	ParFiniteElementSpace fespace(*maker.GetParallelFESpace());
	cout << fespace.GetFE(0) << endl;
	
	int order = 1;
	int nV = pmesh.GetNV();			//number of vertices





	// ======================================
	// FIND DISTANCE FUNCTION BY REINITIALIZING LEVEL SET
	// ======================================
	cout << "FINDING DISTANCE FUNCTION WITH LEVEL SET REINITIALIZATION" << endl;
	
	// Need to use Discontinuous Galerkin (DG) for upwinding (look at MFEM examples 9 and 18)
	DG_FECollection fec_dg(order, pmesh.Dimension(), BasisType::GaussLobatto);
	ParFiniteElementSpace fespace_dg(&pmesh, &fec_dg);
	ParFiniteElementSpace dfespace_dg(&pmesh, &fec_dg, pmesh.Dimension(), Ordering::byNODES); //X1X2X3.....,Y1Y2Y3.....,Z1Z2Z3......
	maker.Make_DG_FESpace_Parallel();

	
	VoxelSolver_DG solver_dg(maker.GetGlobalFESpace(), maker.GetParallelFESpace(), maker.GetParallelFESpace_DG(), maker.GetParallelFESpace_DGdim());
	
	solver_dg.ProjectVals(solver.GetParallelVox());
	solver_dg.CalcLevelSetVel();
		
	ODESolver *ode_solver_dg2 = new ForwardEulerSolver; // WHY IS THIS LINE NECESSARY???
		
	solver_dg.FormMatrices(boundary_dofs);
	// Time Stepping
	t_ode = 0.0;
	dt = 0.1;
	for (int t = 0; t < 10; t++){
		solver_dg.CalcLevelSetVel();
		
		
		solver_dg.UpdateMatricesAndSolve(boundary_dofs, t_ode, dt);
	}
	
	
	// Output Distance to Paraview
	ParGridFunction d(*solver_dg.GetDistFunc());
	//d = d2;
	cout << "PRINTING OUT DistanceFunction" << endl;
	solver.ParaviewSave("DstFun","Dst",&d);

	// Output Advection Velocity to Paraview
	ParGridFunction c(*solver_dg.GetAdvVel());
	//c = c2;
	solver.ParaviewSave("AdvVel","Vel",&c);




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
	solver.ParaviewSave("psi","psi",&psi);
	
	VoxelSolver ConnectSolver(maker.GetGlobalFESpace(),maker.GetParallelFESpace());
		

	for (int ConIter = 0; ConIter < 2; ConIter++){
		ConnectSolver.AssignGlobalValues(0.0);
		ConnectSolver.MapGlobalToLocal(maker.GetGlobalMesh(),maker.GetParallelMesh());
		
		if (ConIter==1){
			psi.Neg();
			psi += 1.0;
		}

		// Dirichlet Boundary conditions
		if (ConIter==0){
			ConnectSolver.NorthDirichletBCs(maker.GetParallelMesh());
		} else{
			ConnectSolver.SouthDirichletBCs(maker.GetParallelMesh());
		}
		ConnectSolver.DetermineConnectivityBCs(psi);
		
		/*
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
		ParGridFunction dbcval(&fespace);
		dbcval = 0.0;
		for (int vi = 0; vi < nV; vi++){
			if ( psi(vi) > 0.6 ){
				dbcval(vi) = 1.0;
			}
		}
		GridFunctionCoefficient dbcCoef(&dbcval);
		
		ConnectSolver.AssignDirichletBCs(dbcCoef, dbc_bdr);
		*/
		
		ConnectSolver.AssignDirichletBCs(*ConnectSolver.GetBCCoef(), *ConnectSolver.GetBCMarker());
		//ConnectSolver.InitMatricesAndTimeDepOpers(ess_tdof_list, psi, psi);
		ConnectSolver.InitMatricesAndTimeDepOpers(*ConnectSolver.GetTDOF(), psi, psi);
		
		t_ode = 0.0;
		dt = 0.05;
		for (int t = 0; t < 100; t++){
			//ConnectSolver.UpdateSystemAndSolve(ess_tdof_list, t_ode, dt);
			//ConnectSolver.AccelerateDiffusion(psi, dbcCoef, dbc_bdr);
			ConnectSolver.UpdateSystemAndSolve(*ConnectSolver.GetTDOF(), t_ode, dt);
			ConnectSolver.AccelerateDiffusion(psi, *ConnectSolver.GetBCCoef(), *ConnectSolver.GetBCMarker());
		}

		// Output Cn to Paraview
		cout << "PRINTING OUT Electrolyte Concentration" << endl;
		if (ConIter==0){
			ConnectSolver.ParaviewSave("Conc_p","Cn_p",ConnectSolver.GetParallelVox());
		} else {
			ConnectSolver.ParaviewSave("Conc_e","Cn_e",ConnectSolver.GetParallelVox());
		}
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
