
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class VoxelSolver {
public:
	VoxelSolver(FiniteElementSpace *fes); //constructor
	void AssignGlobalValues(vector<vector<vector<int>>> data);
	
	GridFunction* GetGlobalVox() {return gVox;}

private:
	GridFunction* gVox;
};
