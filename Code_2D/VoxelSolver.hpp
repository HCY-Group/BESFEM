#ifndef VOXELSOLVER_HPP
#define VOXELSOLVER_HPP

#include "mfem.hpp"
#include "TimeDepOpers.hpp"

//using namespace std;
//using namespace mfem;

class VoxelSolver {
public:
	VoxelSolver(mfem::FiniteElementSpace *fes); //constructor
	VoxelSolver(mfem::FiniteElementSpace *gfes, mfem::ParFiniteElementSpace *fes); //constructor with parallel
	void AssignGlobalValues(double a);
	void AssignGlobalValues(vector<vector<vector<int>>> data);
	void AssignDirichletBCs(mfem::GridFunctionCoefficient &Coef, Array<int> &bdr);
	void ParaviewSave(string FileName, string VariableName, mfem::GridFunction* gf);
	void MapGlobalToLocal(mfem::Mesh* gmesh, mfem::ParMesh* pmesh);
	void InitMatricesAndTimeDepOpers(Array<int> boundary_dofs, mfem::ParGridFunction &Diff, mfem::ParGridFunction &DomPar);
	void UpdateLinearForm(mfem::ParGridFunction gf);
	void UpdateLinearForm_DoubleWellPotential();
	void UpdateSystemAndSolve(Array<int> boundary_dofs, double t_ode, double dt);
	void AccelerateDiffusion(mfem::ParGridFunction &DomPar, mfem::GridFunctionCoefficient &Coef, Array<int> &bdr);
	void NorthDirichletBCs(mfem::Mesh *mesh);
	void SouthDirichletBCs(mfem::Mesh *mesh);
	void DetermineConnectivityBCs(mfem::ParGridFunction &DomPar);
	void MultiplyVox(mfem::GridFunction &gf); //Elementwise Multiplication
	
	mfem::GridFunction* GetGlobalVox() {return gVox;}
	mfem::ParGridFunction* GetParallelVox() {return Vox;}
	
	mfem::GridFunctionCoefficient* GetBCCoef() {return dbcCoef;}
	Array<int>* GetBCMarker() {return dbc_bdr;}
	Array<int>* GetTDOF() {return ess_tdof_list;}

protected:
	mfem::ODESolver* ode_solver;
	
	Array<int> *ess_tdof_list = nullptr;	
	Array<int> *dbc_bdr = nullptr;
	mfem::ParGridFunction *dbcval = nullptr;
	mfem::GridFunctionCoefficient *dbcCoef = nullptr;
	
private:
	mfem::GridFunction* gVox;
	mfem::ParGridFunction* Vox;
	
	mfem::ParBilinearForm *K;
	mfem::ParLinearForm *Fct;
	
	ConductionOperator* oper;
};

#endif
