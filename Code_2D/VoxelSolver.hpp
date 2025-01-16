#ifndef VOXELSOLVER_HPP
#define VOXELSOLVER_HPP

#include "mfem.hpp"
#include "TimeDepOpers.hpp"

using namespace std;
using namespace mfem;

class VoxelSolver {
public:
	VoxelSolver(FiniteElementSpace *fes); //constructor
	VoxelSolver(FiniteElementSpace *gfes, ParFiniteElementSpace *fes); //constructor with parallel
	void AssignGlobalValues(double a);
	void AssignGlobalValues(vector<vector<vector<int>>> data);
	void AssignDirichletBCs(GridFunctionCoefficient &Coef, Array<int> &bdr);
	void ParaviewSave(string FileName, string VariableName, GridFunction* gf);
	void MapGlobalToLocal(Mesh* gmesh, ParMesh* pmesh);
	//void InitStiffMatrix(Array<int> boundary_dofs, ParGridFunction Diff);
	//void InitTimeDepOper(ParGridFunction DomPar);
	void InitMatricesAndTimeDepOpers(Array<int> boundary_dofs, ParGridFunction &Diff, ParGridFunction &DomPar);
	void UpdateLinearForm(ParGridFunction gf);
	void UpdateLinearForm_DoubleWellPotential();
	void UpdateSystemAndSolve(Array<int> boundary_dofs, double t_ode, double dt);
	void AccelerateDiffusion(ParGridFunction &DomPar, GridFunctionCoefficient &Coef, Array<int> &bdr);
	void NorthDirichletBCs(Mesh *mesh);
	void SouthDirichletBCs(Mesh *mesh);
	void DetermineConnectivityBCs(ParGridFunction &DomPar);
	
	GridFunction* GetGlobalVox() {return gVox;}
	ParGridFunction* GetParallelVox() {return Vox;}
	
	GridFunctionCoefficient* GetBCCoef() {return dbcCoef;}
	Array<int>* GetBCMarker() {return dbc_bdr;}
	Array<int>* GetTDOF() {return ess_tdof_list;}

protected:
	ODESolver* ode_solver;
	
	Array<int> *ess_tdof_list = nullptr;	
	Array<int> *dbc_bdr = nullptr;
	ParGridFunction *dbcval = nullptr;
	GridFunctionCoefficient *dbcCoef = nullptr;
	
private:
	GridFunction* gVox;
	ParGridFunction* Vox;
	
	ParBilinearForm *K;
	ParLinearForm *Fct;
	//HypreParMatrix Kmat;
	//HypreParVector Fcb;
	
	ConductionOperator* oper;
};

#endif
