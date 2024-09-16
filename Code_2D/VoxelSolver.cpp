
#include "VoxelSolver.hpp"

using namespace mfem;
using namespace std;

VoxelSolver::VoxelSolver(FiniteElementSpace *fes){
	gVox = new GridFunction(fes);
	cout << "GVOX SIZE" << gVox->Size() << endl;
}
	
void VoxelSolver::AssignGlobalValues(vector<vector<vector<int>>> data) {
	// Create global FE space for Voxel Data
	/*
	int order = 1;
	H1_FECollection gFec(order, gmesh.Dimension());
	gFespace = new FiniteElementSpace(&gmesh, &gFec);
	*/

	// global grid function for voxel data
	cout << "Defining Voxel GridFunction" << endl;
	//gVox = new GridFunction(gFespace);
	int nz = data.size();
	int ny = data[0].size();
	int nx = data[0][0].size();
	Vector tmp(gVox->Size());
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				int idx = i + nx*j + nx*ny*k;
				//cout << idx << endl;
				//cout << "idx = " << i+ny*j << endl;
				//cout << "gVox[idx] = " << gVox[idx] << endl;
				//cout << "data[i][j][0] = " << data[i][j][0] << endl;
				
				//gVox[idx] = data[k][j][i];
				tmp[idx] = data[k][j][i];
			}
		}
	}
	*gVox = tmp;

}
