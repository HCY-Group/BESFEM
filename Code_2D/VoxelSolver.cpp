
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
	//Vector tmp(gVox->Size());
	GridFunction tmp(*gVox); // Copy constructor (using to "dereference")
	tmp = *gVox;
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
	//cout << "gVox fec: " << gVox->OwnFEC() << endl;
	
	//FiniteElementSpace* fes = gVox->FESpace();
}

void VoxelSolver::ParaviewSave(string FileName, string VariableName, GridFunction* gf) {
	
	// weirdness to get ParaView save to work
	// Can't just return the grid function
	// We have to create the FiniteElementSpace again from scratch, starting with the FECollection
	
	int order = 1;

	GridFunction gf3(*gf);
	FiniteElementSpace* fes = gf->FESpace();
	Mesh* mesh = fes->GetMesh();
	cout << "MESH " << mesh << endl;
	cout << "GF " << gf << endl;
	//const FiniteElementCollection* gFec = fes->FEColl();
	//cout << "gVox fec: " << gf->OwnFEC() << endl;
	//cout << "fes fec: " << fes->FEColl() << endl;
	
	H1_FECollection fec2(order, mesh->Dimension());
	FiniteElementSpace fes2(mesh, &fec2);
	GridFunction gf2(&fes2);
	
	gf2 = gf3;
	
	ParaViewDataCollection *pd = NULL;
	cout << "before" << endl;
	//pd = new ParaViewDataCollection("gVoxelData", &gmesh2);
	
	
	//gf->MakeOwner(fes->FEColl());
	
	pd = new ParaViewDataCollection(FileName, mesh);
	//pd->RegisterField("gVox", &gVox);
	//pd->RegisterField(VariableName, gf);
	pd->RegisterField(VariableName, &gf2);
	//pd->RegisterField(VariableName, gVox);
	//pd->RegisterField("gVox", maker.GetGlobalVox());
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	cout << "after" << endl;
	pd->Save();
	cout << "after" << endl;
	delete pd;
}
