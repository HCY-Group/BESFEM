
#include "mfem.hpp"
#include "TimeDepOpers.hpp"

using namespace std;
using namespace mfem;

class VoxelSolver {
public:
	VoxelSolver(FiniteElementSpace *fes); //constructor
	VoxelSolver(FiniteElementSpace *gfes, ParFiniteElementSpace *fes); //constructor with parallel
	void AssignGlobalValues(vector<vector<vector<int>>> data);
	void ParaviewSave(string FileName, string VariableName, GridFunction* gf);
	void MapGlobalToLocal(Mesh* gmesh, ParMesh* pmesh);
	//void InitStiffMatrix(Array<int> boundary_dofs, ParGridFunction Diff);
	//void InitTimeDepOper(ParGridFunction DomPar);
	void InitMatricesAndTimeDepOpers(Array<int> boundary_dofs, ParGridFunction Diff, ParGridFunction DomPar);
	void UpdateLinearForm(ParGridFunction gf);
	void UpdateLinearForm_DoubleWellPotential();
	void UpdateSystemAndSolve(Array<int> boundary_dofs, double t_ode, double dt);
	
	GridFunction* GetGlobalVox() {return gVox;}
	ParGridFunction* GetParallelVox() {return Vox;}
	
private:
	GridFunction* gVox;
	ParGridFunction* Vox;
	
	ParBilinearForm *K;
	ParLinearForm *Fct;
	//HypreParMatrix Kmat;
	//HypreParVector Fcb;
	
	ConductionOperator* oper;
	ODESolver* ode_solver;
};
