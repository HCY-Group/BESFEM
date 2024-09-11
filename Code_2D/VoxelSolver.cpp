
#include "VoxelSolver.hpp"

using namespace mfem;
using namespace std;

VoxelSolver::VoxelSolver(FiniteElementSpace *fes){
	gVox = new GridFunction(fes);
	cout << "GVOX SIZE" << gVox->Size() << endl;
}
	
