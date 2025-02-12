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
	//TODO: The code works with serial 2d, parallel 2d, and serial 3d, but not parallel 3d
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







	// ======================================
	// FIND DISTANCE FUNCTION BY REINITIALIZING LEVEL SET
	// ======================================
	cout << "FINDING DISTANCE FUNCTION WITH LEVEL SET REINITIALIZATION" << endl;
	
	// Need to use Discontinuous Galerkin (DG) for upwinding (look at MFEM examples 9 and 18)
	int order = 1;
	//DG_FECollection fec_dg(order, maker.GetParallelMesh()->Dimension(), BasisType::GaussLobatto);
	//ParFiniteElementSpace fespace_dg(maker.GetParallelMesh(), &fec_dg);
	//ParFiniteElementSpace dfespace_dg(maker.GetParallelMesh(), &fec_dg, maker.GetParallelMesh()->Dimension(), Ordering::byNODES); //X1X2X3.....,Y1Y2Y3.....,Z1Z2Z3......
	maker.Make_DG_FESpace_Parallel();

	
	VoxelSolver_DG solver_dg(maker.GetGlobalFESpace(), maker.GetParallelFESpace(), maker.GetParallelFESpace_DG(), maker.GetParallelFESpace_DGdim());
	
	solver_dg.ProjectVals(solver.GetParallelVox());
	solver_dg.CalcLevelSetVel();
		
	ODESolver *ode_solver_dg2 = new ForwardEulerSolver; //TODO: WHY IS THIS LINE NECESSARY???
		
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
	ParGridFunction psi(maker.GetParallelFESpace());
	psi.ProjectGridFunction(d);
	psi -= 0.5; // Center about 0
	for (int vi = 0; vi < maker.GetParallelMesh()->GetNV(); vi++){
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
		
		
		ConnectSolver.AssignDirichletBCs(*ConnectSolver.GetBCCoef(), *ConnectSolver.GetBCMarker());
		ConnectSolver.InitMatricesAndTimeDepOpers(*ConnectSolver.GetTDOF(), psi, psi);
		
		t_ode = 0.0;
		dt = 0.05;
		//dt = 0.1;
		for (int t = 0; t < 300; t++){
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




	return 0;

}
