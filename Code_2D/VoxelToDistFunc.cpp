#include "mfem.hpp"
//#include "myHeatDistanceSolver.hpp"
#include "dist_solver.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include "TimeDepOpers.hpp"
#include "readtiff.h"
#include "MeshMaker.hpp"
#include "VoxelSolver.hpp"
#include "VoxelSolver_DG.hpp"
#include "../code/Initialize_Geometry.hpp"
#include "../code/Constants.hpp"

using namespace std;
using namespace mfem;



void ComputeScalarDistance(ParGridFunction &source, ParGridFunction &distance);
void MyComputeDistance(ParGridFunction &psi, ParGridFunction &distance);

int main(int argc, char *argv[])
{
	// Initialize MPI and HYPRE
	Mpi::Init(argc, argv);
	Hypre::Init();

	/*	
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
	
	*/
	
	
	






	// ======================================
	// LOCAL (PARALLEL) MESH AND DATA
	// ======================================
	

	// Test Initialize_Geometry class
    	Initialize_Geometry geometry;
    	geometry.InitializeMesh(Constants::mesh_file, MPI_COMM_WORLD, Constants::order);
    	geometry.SetupBoundaryConditions();

	//VoxelSolver solver(maker.GetGlobalFESpace());
	//VoxelSolver solver(maker.GetGlobalFESpace(), maker.GetParallelFESpace());
	VoxelSolver solver(&*geometry.globalfespace, &*geometry.parfespace);
	//solver.AssignGlobalValues(tiffdata);
	solver.AssignGlobalValues(geometry.tiffData);
	//solver.ParaviewSave("gVoxelData","gVox",solver.GetGlobalVox());

	//solver.MapGlobalToLocal(maker.GetGlobalMesh(),maker.GetParallelMesh());
	solver.MapGlobalToLocal(&*geometry.globalMesh,&*geometry.parallelMesh);
	solver.ParaviewSave("pVoxelData","Vox",solver.GetParallelVox());
	

	ODESolver *ode_solver_dg2 = new ForwardEulerSolver; // WHY IS THIS LINE NECESSARY???
	
	// ======================================
	// SMOOTH THE DATA USING ALLEN-CAHN
	// ======================================
	cout << "SMOOTHING WITH ALLEN-CAHN" << endl;
	
	// Natural (Neumann) boundary conditions
	Array<int> boundary_dofs;
	
	// Domain Parameter
	//ParGridFunction ones(maker.GetParallelFESpace());
	ParGridFunction ones(&*geometry.parfespace);
	ones = 1.0;
	
	// Diffusion Coefficient
	//ParGridFunction Mob(&fespace);
	//ParGridFunction Mob(maker.GetParallelFESpace());
	ParGridFunction Mob(&*geometry.parfespace);
	Mob = 0.64;
	

	solver.InitMatricesAndTimeDepOpers(boundary_dofs, Mob, ones);
	
	// time step
	double t_ode = 0.0;
	double dt = 0.05;
	for (int t = 0; t < 20; t++){

		solver.UpdateLinearForm_DoubleWellPotential();
		
		solver.UpdateSystemAndSolve(boundary_dofs, t_ode, dt);
		
		double lVoxMax, lVoxMin, gVoxMax, gVoxMin;
		lVoxMax = solver.GetParallelVox()->Max();
		lVoxMin = solver.GetParallelVox()->Min();
		MPI_Reduce(&lVoxMax, &gVoxMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lVoxMin, &gVoxMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		if (mfem::Mpi::WorldRank()==0){
			cout << "Max Vox: " << gVoxMax << " Min Vox: " << gVoxMin << endl;
		}
	}
	//oper.~ConductionOperator();
	
	// Output Vox to Paraview
	solver.ParaviewSave("SmoothVox","Vox",solver.GetParallelVox());








	// ======================================
	// FIND DISTANCE FUNCTION using HeatDistancing in MFEM
	// See the Distancing miniapp in MFEM
	// This method is based on the article
	// K. Crane, C. Weischedel, M. Weischedel
	// Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow
	// ACM Transactions on Graphics, Vol. 32, No. 5, October, 2013
	// ======================================

	/*
	GridFunctionCoefficient ls_coeff(solver.GetParallelVox());
	const real_t dx = geometry.GetParallelMesh()->GetElementSize(0);
	real_t t_param = 1.0;
	//DistanceSolver *dist_solver = NULL;
	//auto ds = new HeatDistanceSolver(t_param*dx*dx);
	HeatDistanceSolver ds(t_param*dx*dx);
	//DistanceSolver *ds2 = new HeatDistanceSolver(t_param*dx*dx);
	ParGridFunction distance_s(&*geometry.GetParFiniteElementSpace());
	//ds.ComputeScalarDistance(ls_coeff, distance_s);  //<--- This line is giving 'vtable' errors...
	*/

	/*
	// Define initial conditions
	ParGridFunction u(*solver.GetParallelVox());
	// Bump-like source term at boundary
	for (int i = 0; i < u.Size(); i++){
		real_t x = u(i);
		u(i) = 4.0 * (1.0-x) * x;
	}

	ParGridFunction ds(u);
	ds = 0.0;
	ComputeScalarDistance(u, ds);
   
	cout << "MAX ds: " << ds.Max() << " MIN ds: " << ds.Min() << endl;
	solver.ParaviewSave("HeatDst","Dst",&ds);
	*/

	// Find distance using my own...
	cout << "My own distancing function" << endl;
	
	/*
	ParGridFunction ds(*solver.GetParallelVox());
	ds = 0.0;
	MyComputeDistance(*solver.GetParallelVox(), ds);
	*/
	
	//Normalize psi
	//*solver.GetParallelVox() -= solver.GetParallelVox()->Min();
	//*solver.GetParallelVox() += 1e-7;
	/*
	for (int i = 0; i < solver.GetParallelVox()->Size(); i++) {
		if ( (*solver.GetParallelVox())(i) < 1e-5 ) {(*solver.GetParallelVox())(i) = 1e-5;}
		if ( (*solver.GetParallelVox())(i) > 1.00 ) {(*solver.GetParallelVox())(i) = 1.00;}
	}
	*/
	VoxelSolver solver_dist(&*geometry.globalfespace, &*geometry.parfespace);
	solver_dist.AssignGlobalValues(0.0);
	//solver_dist.AssignGlobalValues(1.0);
	solver_dist.MapGlobalToLocal(&*geometry.globalMesh,&*geometry.parallelMesh);

	//solver_dist.InitMatricesAndTimeDepOpers(boundary_dofs, *solver.GetParallelVox(), *solver.GetParallelVox());
	solver_dist.InitMatricesAndTimeDepOpers(boundary_dofs, *solver.GetParallelVox(), ones);
	
	// time step
	t_ode = 0.0;
	dt = 1e-5;
	dt = 0.05;
	for (int t = 0; t < 100; t++){

		solver_dist.UpdateLinearForm_SBMDirichlet(*solver.GetParallelVox());
		
		solver_dist.UpdateSystemAndSolve(boundary_dofs, t_ode, dt);
		
		double lVoxMax, lVoxMin, gVoxMax, gVoxMin;
		lVoxMax = solver_dist.GetParallelVox()->Max();
		lVoxMin = solver_dist.GetParallelVox()->Min();
		MPI_Reduce(&lVoxMax, &gVoxMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lVoxMin, &gVoxMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		if (mfem::Mpi::WorldRank()==0){
			cout << "Max Vox: " << gVoxMax << " Min Vox: " << gVoxMin << endl;
		}
	}
	//oper.~ConductionOperator();
	
	// Output Vox to Paraview
	solver_dist.ParaviewSave("PoissonVox","C",solver_dist.GetParallelVox());
	

	












	
	// ======================================
	// FIND DISTANCE FUNCTION BY REINITIALIZING LEVEL SET
	// ======================================
	cout << "FINDING DISTANCE FUNCTION WITH LEVEL SET REINITIALIZATION" << endl;
	
	// Need to use Discontinuous Galerkin (DG) for upwinding (look at MFEM examples 9 and 18)

	
	//VoxelSolver_DG solver_dg(maker.GetGlobalFESpace(), maker.GetParallelFESpace(), maker.GetParallelFESpace_DG(), maker.GetParallelFESpace_DGdim());
	VoxelSolver_DG solver_dg(&*geometry.globalfespace, &*geometry.parfespace, &*geometry.parfespace_dg, &*geometry.pardimfespace_dg);
	
	solver_dg.ProjectVals(solver.GetParallelVox());
	solver_dg.CalcLevelSetVel();
		
	solver_dg.FormMatrices(boundary_dofs);
	// Time Stepping
	t_ode = 0.0;
	dt = 0.01;
	//for (int t = 0; t < 500; t++){
	for (int t = 0; t < 200; t++){
		solver_dg.CalcLevelSetVel();
		
		
		solver_dg.UpdateMatricesAndSolve(boundary_dofs, t_ode, dt);
		
		double lVoxMax, lVoxMin, gVoxMax, gVoxMin;
		lVoxMax = solver_dg.GetDistFunc()->Max();
		lVoxMin = solver_dg.GetDistFunc()->Min();
		MPI_Reduce(&lVoxMax, &gVoxMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lVoxMin, &gVoxMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		if (mfem::Mpi::WorldRank()==0){
			cout << "Max d: " << gVoxMax << " Min d: " << gVoxMin << endl;
		}
	}
	
	
	// Output Distance to Paraview
	ParGridFunction d(*solver_dg.GetDistFunc());
	d -= 0.5; // Center about 0
	//d = d2;
	cout << "PRINTING OUT DistanceFunction" << endl;
	solver.ParaviewSave("DstFun","Dst",&d);

	// Output Advection Velocity to Paraview
	ParGridFunction c(*solver_dg.GetAdvVel());
	//c = c2;
	//solver.ParaviewSave("AdvVel","Vel",&c);


	
	
	// ======================================
	// FIND CONNECTIVITY
	// ======================================

	// Use d to define psi
	//ParGridFunction psi(maker.GetParallelFESpace());
	ParGridFunction psi(&*geometry.parfespace);
	psi.ProjectGridFunction(d);
	//psi -= 0.5; // Center about 0
	//for (int vi = 0; vi < maker.GetParallelMesh()->GetNV(); vi++){
	for (int vi = 0; vi < geometry.parallelMesh->GetNV(); vi++){
		psi(vi) = 0.5*( tanh(psi(vi)) + 1.0 );
	}

	// Output psi to Paraview
	cout << "PRINTING OUT Electrolyte Concentration" << endl;
	solver.ParaviewSave("psi","psi",&psi);
	
		

	for (int ConIter = 0; ConIter < 2; ConIter++){
		//VoxelSolver ConnectSolver(maker.GetGlobalFESpace(),maker.GetParallelFESpace());
		//VoxelSolver_DG solver_dg2(maker.GetGlobalFESpace(), maker.GetParallelFESpace(), maker.GetParallelFESpace_DG(), maker.GetParallelFESpace_DGdim());
		VoxelSolver ConnectSolver(&*geometry.globalfespace,&*geometry.parfespace);
		VoxelSolver_DG solver_dg2(&*geometry.globalfespace, &*geometry.parfespace, &*geometry.parfespace_dg, &*geometry.pardimfespace_dg);
		
		ConnectSolver.AssignGlobalValues(0.0);
		//ConnectSolver.MapGlobalToLocal(maker.GetGlobalMesh(),maker.GetParallelMesh());
		ConnectSolver.MapGlobalToLocal(&*geometry.globalMesh,&*geometry.parallelMesh);
		
		if (ConIter==1){
			psi.Neg();
			psi += 1.0;
		}

		// Dirichlet Boundary conditions
		if (ConIter==0){
			//ConnectSolver.NorthDirichletBCs(maker.GetParallelMesh());
			ConnectSolver.NorthDirichletBCs(&*geometry.parallelMesh);
		} else{
			//ConnectSolver.SouthDirichletBCs(maker.GetParallelMesh());
			ConnectSolver.SouthDirichletBCs(&*geometry.parallelMesh);
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
		
			double lVoxMax, lVoxMin, gVoxMax, gVoxMin;
			lVoxMax = ConnectSolver.GetParallelVox()->Max();
			lVoxMin = ConnectSolver.GetParallelVox()->Min();
			MPI_Reduce(&lVoxMax, &gVoxMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&lVoxMin, &gVoxMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
			if (mfem::Mpi::WorldRank()==0){
				cout << "Max Vox: " << gVoxMax << " Min Vox: " << gVoxMin << endl;
			}
		}
		
		//Multiply value by psi so the boundaries are closer to each other
		ConnectSolver.MultiplyVox(psi);
		
		// Output Cn to Paraview
		cout << "PRINTING OUT Electrolyte Concentration" << endl;
		if (ConIter==0){
			ConnectSolver.ParaviewSave("Conc_p","Cn_p",ConnectSolver.GetParallelVox());
		} else {
			ConnectSolver.ParaviewSave("Conc_e","Cn_e",ConnectSolver.GetParallelVox());
		}

		// ======================================
		// FIND DISTANCE FUNCTION BY REINITIALIZING LEVEL SET
		// ======================================
		cout << "FINDING DISTANCE FUNCTION WITH LEVEL SET REINITIALIZATION" << endl;
		
		solver_dg2.ProjectVals(ConnectSolver.GetParallelVox());
		solver_dg2.CalcLevelSetVel();
			
		//ODESolver *ode_solver_dg2 = new ForwardEulerSolver; //TODO: WHY IS THIS LINE NECESSARY???
			
		solver_dg2.FormMatrices(boundary_dofs);
		// Time Stepping
		t_ode = 0.0;
		dt = 0.01;
		//for (int t = 0; t < 500; t++){
		for (int t = 0; t < 200; t++){
			solver_dg2.CalcLevelSetVel();
			
			
			solver_dg2.UpdateMatricesAndSolve(boundary_dofs, t_ode, dt);
		
			double lVoxMax, lVoxMin, gVoxMax, gVoxMin;
			lVoxMax = solver_dg2.GetDistFunc()->Max();
			lVoxMin = solver_dg2.GetDistFunc()->Min();
			MPI_Reduce(&lVoxMax, &gVoxMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&lVoxMin, &gVoxMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
			if (mfem::Mpi::WorldRank()==0){
				cout << "Max d: " << gVoxMax << " Min d: " << gVoxMin << endl;
			}
		}
		
		
		// Output Distance to Paraview
		ParGridFunction d2(*solver_dg2.GetDistFunc());
		//d = d2;
	
		// Project d2 onto H1 FESpace to use in BESFEM sim
		//ParGridFunction psi2(maker.GetParallelFESpace());
		ParGridFunction psi2(&*geometry.parfespace);
		psi2.ProjectGridFunction(d2);
		psi2 -= 0.5; // Center about 0
		
		cout << "PRINTING OUT DistanceFunction" << endl;
		if (ConIter==0){
			//solver.ParaviewSave("DstFun_p","Dst_p",&d2);
			solver.ParaviewSave("DstFun_p","Dst_p",&psi2);
			psi2.SaveAsOne("dsF_p.txt");

			// Output Advection Velocity to Paraview
			ParGridFunction c(*solver_dg2.GetAdvVel());
			//c = c2;
			//solver.ParaviewSave("AdvVel_p","Vel_p",&c);
		} else {
			//solver.ParaviewSave("DstFun_e","Dst_e",&d2);
			solver.ParaviewSave("DstFun_e","Dst_e",&psi2);
			psi2.SaveAsOne("dsF_e.txt");

			// Output Advection Velocity to Paraview
			ParGridFunction c(*solver_dg.GetAdvVel());
			//c = c2;
			//solver.ParaviewSave("AdvVel_e","Vel_e",&c);
		}

	}

	// Output mesh to use in BESFEM sims
	//Mesh mesh2(*maker.GetGlobalMesh());
	Mesh mesh2(*geometry.globalMesh);
	mesh2.Save("VoxMesh.mesh");
	


	return 0;

}



/*
class NormalizedGradCoefficient : public VectorCoefficient
{
private:
   const GridFunction &u;

public:
   NormalizedGradCoefficient(const GridFunction &u_gf, int dim)
      : VectorCoefficient(dim), u(u_gf) { }

   using VectorCoefficient::Eval;

   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
   {
      T.SetIntPoint(&ip);

      u.GetGradient(T, V);
      const double norm = V.Norml2() + 1e-12;
      V /= -norm;
   }
};
*/



void ComputeScalarDistance(ParGridFunction &source, ParGridFunction &distance){
   
   int diffuse_iter = 1;
   real_t parameter_t = 1;
   int cg_print_lvl = 0;
   int amg_print_lvl = 0;

   ParFiniteElementSpace &pfes = *distance.ParFESpace();
   ParMesh &pmesh = *pfes.GetParMesh();

   // Solver.
   CGSolver cg(MPI_COMM_WORLD);
   cg.SetRelTol(1e-12);
   cg.SetMaxIter(100);
   cg.SetPrintLevel(cg_print_lvl);
   OperatorPtr A;
   Vector B, X;

   // Step 1 - diffuse.
   ParGridFunction diffused_source(&pfes);
   for (int i = 0; i < diffuse_iter; i++)
   {
      // Set up RHS.
      ParLinearForm b(&pfes);
      GridFunctionCoefficient src_coeff(&source);
      b.AddDomainIntegrator(new DomainLFIntegrator(src_coeff));
      b.Assemble();

      // Diffusion and mass terms in the LHS.
      ParBilinearForm a_d(&pfes);
      a_d.AddDomainIntegrator(new MassIntegrator);
      ConstantCoefficient t_coeff(parameter_t);
      a_d.AddDomainIntegrator(new DiffusionIntegrator(t_coeff));
      a_d.Assemble();

      // Solve with Dirichlet BC.
      Array<int> ess_tdof_list;
      if (pmesh.bdr_attributes.Size())
      {
         Array<int> ess_bdr(pmesh.bdr_attributes.Max());
         ess_bdr = 1;
         pfes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
      }
      ParGridFunction u_dirichlet(&pfes);
      u_dirichlet = 0.0;
      a_d.FormLinearSystem(ess_tdof_list, u_dirichlet, b, A, X, B);
      auto *prec = new HypreBoomerAMG;
      prec->SetPrintLevel(amg_print_lvl);
      cg.SetPreconditioner(*prec);
      cg.SetOperator(*A);
      cg.Mult(B, X);
      a_d.RecoverFEMSolution(X, b, u_dirichlet);
      delete prec;

	cout << "HERE, end of Dirichlet" << endl;
	cout << "MAX u: " << u_dirichlet.Max() << " MIN u: " << u_dirichlet.Min() << endl;

      // Diffusion and mass terms in the LHS.
      ParBilinearForm a_n(&pfes);
      a_n.AddDomainIntegrator(new MassIntegrator);
      a_n.AddDomainIntegrator(new DiffusionIntegrator(t_coeff));
      a_n.Assemble();

      // Solve with Neumann BC.
      ParGridFunction u_neumann(&pfes);
      ess_tdof_list.DeleteAll();
      a_n.FormLinearSystem(ess_tdof_list, u_neumann, b, A, X, B);
      auto *prec2 = new HypreBoomerAMG;
      prec2->SetPrintLevel(amg_print_lvl);
      cg.SetPreconditioner(*prec2);
      cg.SetOperator(*A);
      cg.Mult(B, X);
      a_n.RecoverFEMSolution(X, b, u_neumann);
      delete prec2;

	cout << "HERE, end of Neumann" << endl;
	cout << "MAX u: " << u_neumann.Max() << " MIN u: " << u_neumann.Min() << endl;

      for (int ii = 0; ii < diffused_source.Size(); ii++)
      {
         // This assumes that the magnitudes of the two solutions are somewhat
         // similar; otherwise one of the solutions would dominate and the BC
         // won't look correct. To avoid this, it's good to have the source
         // away from the boundary (i.e. have more resolution).
         diffused_source(ii) = 0.5 * (u_neumann(ii) + u_dirichlet(ii));
         diffused_source(ii) = u_neumann(ii);
      }
      source = diffused_source;
   }
   

	ParaViewDataCollection *pd = NULL;
	pd = new ParaViewDataCollection("dist_out1", diffused_source.FESpace()->GetMesh());
	pd->RegisterField("diffsource", &diffused_source);
	pd->SetLevelsOfDetail(1);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;


   // Step 2 - solve for the distance using the normalized gradient.
   {
      // RHS - normalized gradient.
      ParLinearForm b2(&pfes);
      NormalizedGradCoefficient grad_u(diffused_source, pmesh.Dimension());
      b2.AddDomainIntegrator(new DomainLFGradIntegrator(grad_u));
      b2.Assemble();

      // LHS - diffusion.
      ParBilinearForm a2(&pfes);
      a2.AddDomainIntegrator(new DiffusionIntegrator);
      a2.Assemble();

      // No BC.
      Array<int> no_ess_tdofs;

      a2.FormLinearSystem(no_ess_tdofs, distance, b2, A, X, B);

      auto *prec = new HypreBoomerAMG;
      prec->SetPrintLevel(amg_print_lvl);
      cg.SetPreconditioner(*prec);
      cg.SetOperator(*A);
      cg.Mult(B, X);
      a2.RecoverFEMSolution(X, b2, distance);
      delete prec;
   }

	cout << "HERE, end of Distance" << endl;
	cout << "MAX u: " << distance.Max() << " MIN u: " << distance.Min() << endl;

	//ParaViewDataCollection *pd = NULL;
	pd = new ParaViewDataCollection("dist_out2", distance.FESpace()->GetMesh());
	pd->RegisterField("dist", &distance);
	pd->SetLevelsOfDetail(1);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;
   
   // Shift the distance values to have minimum at zero.
   double d_min_loc = distance.Min();
   double d_min_glob;
   MPI_Allreduce(&d_min_loc, &d_min_glob, 1, MPI_DOUBLE,
                 MPI_MIN, pfes.GetComm());
   distance -= d_min_glob;
}



void MyComputeDistance(ParGridFunction &psi, ParGridFunction &distance){
   cout << "Min psi: " << psi.Min() << ", Max psi: " << psi.Max() << endl;
   ParFiniteElementSpace &pfes = *distance.ParFESpace();
   ParMesh &pmesh = *pfes.GetParMesh();

   // Solve Lap(u) = -1, with Dirichlet BC u = 0
   // A. Belyaev et al: "On Variational and PDE-based Distance Function Approximations", Section 7
   // SBM formulation: div(psi grad(u)) - grad(psi)*grad(u) = -psi
   // SBM formulation: -div(psi grad(u)) + grad(psi) dot grad(u) = psi
   
   for (int i = 0; i < psi.Size(); i++){
      if (psi(i) < 1e-7) {psi(i) = 1e-7;}
   }
	ParaViewDataCollection *pd = NULL;
	pd = new ParaViewDataCollection("MyDistTest_psi", psi.FESpace()->GetMesh());
	pd->RegisterField("u", &psi);
	pd->SetLevelsOfDetail(1);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;
   // Set RHS
   ParLinearForm b(&pfes);
   GridFunctionCoefficient psi_coef(&psi);
   b.AddDomainIntegrator(new DomainLFIntegrator(psi_coef));
   b.Assemble();
	cout << "max b: " << b.Max() << " min b: " << b.Min() << endl;

   // Set diffusion operator
   ParBilinearForm a(&pfes);
   a.AddDomainIntegrator(new DiffusionIntegrator(psi_coef));
   //ParMixedBilinearForm a(&pfes);
   //a.AddDomainIntegrator(new MixedGradGradIntegrator(psi_coef));
   GradientGridFunctionCoefficient psigrad_coef(&psi);
   a.AddDomainIntegrator(new ConvectionIntegrator(psigrad_coef));
   a.Assemble();

   //OperatorPtr A;
   HypreParMatrix A;
   HypreParVector B, X;
   Array<int> ess_tdof_list;
   ParGridFunction x(&pfes);
   //x = psi;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

/*
   HypreSmoother prec;
   prec.SetType(HypreSmoother::Jacobi);
   //HypreSolver *amg = new HypreBoomerAMG;
   //CGSolver cg(MPI_COMM_WORLD);
   //GMRESSolver cg(MPI_COMM_WORLD);
   FGMRESSolver cg(MPI_COMM_WORLD);
   //BiCGSTABSolver cg(MPI_COMM_WORLD);
   cg.SetRelTol(1e-12);
   cg.SetMaxIter(3000);
   //cg.SetPreconditioner(prec);
   //cg.SetPreconditioner(*amg);
   cg.SetOperator(A);
   cg.SetPrintLevel(3);
   cg.Mult(B, X);
   a.RecoverFEMSolution(X, b, x);

   //x -= x.Min();
   //x *= psi;
*/
	// solve for psi dC/dt
	HypreParVector dCdt(X);
	A.Mult(X,dCdt);
	cout << "here" << endl;
	dCdt += B;
	cout << "here" << endl;
	a.RecoverFEMSolution(dCdt, b, x);
	

	cout << "My x max: " << x.Max() << " My x min: " << x.Min() << endl;

	//ParaViewDataCollection *pd = NULL;
	pd = new ParaViewDataCollection("MyDistTest", x.FESpace()->GetMesh());
	//pd->RegisterField("u", &x);
	pd->RegisterField("dCdt", &x);
	pd->SetLevelsOfDetail(1);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;





}













