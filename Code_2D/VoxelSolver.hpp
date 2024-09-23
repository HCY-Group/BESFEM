
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class VoxelSolver {
public:
	VoxelSolver(FiniteElementSpace *fes); //constructor
	void AssignGlobalValues(vector<vector<vector<int>>> data);
	void ParaviewSave(string FileName, string VariableName, GridFunction* gf);
	
	GridFunction* GetGlobalVox() {return gVox;}

private:
	GridFunction* gVox;
};
