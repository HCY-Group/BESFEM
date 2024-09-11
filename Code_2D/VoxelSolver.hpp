
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class VoxelSolver {
public:
	VoxelSolver(FiniteElementSpace *fes); //constructor
	
	GridFunction* GetGlobal() {return gVox;}

private:
	GridFunction* gVox;
};
