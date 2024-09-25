
#include "MeshMaker.hpp"

using namespace mfem;
using namespace std;

MeshMaker::MeshMaker(std::vector<std::vector<std::vector<int>>> tiffdata){
	data = tiffdata;
}

void MeshMaker::MakeGlobalMesh() {
	int nz = data.size();
	int ny = data[0].size();
	int nx = data[0][0].size();
	double sx = nx;  //make dx = 1
	double sy = ny;  //make dy = 1
	double sz = nz;  //make dz = 1
	bool generate_edges = false;
	bool sfc_ordering = false;
	//Mesh gmesh = new Mesh;
	if (nz == 1) {
		gmesh = Mesh::MakeCartesian2D(nx-1, ny-1, Element::QUADRILATERAL, generate_edges, sx, sy, sfc_ordering);
	} else {
		gmesh = Mesh::MakeCartesian3D(nx-1, ny-1, nz-1, Element::HEXAHEDRON, sx, sy, sz, sfc_ordering);
	}
	gmesh.EnsureNCMesh(true);

	cout << "Mesh from MeshMaker " << &gmesh << endl;
	/*
	// Create global FE space for Voxel Data
	int order = 1;
	H1_FECollection gFec(order, gmesh.Dimension());
	FiniteElementSpace gFespace(&gmesh, &gFec);

	// global grid function for voxel data
	cout << "Defining Voxel GridFunction" << endl;
	//GridFunction gVox(&gFespace);
	gVox.SetSpace(&gFespace);
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				int idx = i + nx*j + nx*ny*k;
				//cout << "idx = " << i+ny*j << endl;
				//cout << "gVox[idx] = " << gVox[idx] << endl;
				//cout << "data[i][j][0] = " << data[i][j][0] << endl;
				gVox[idx] = data[k][j][i];
			}
		}
	}
	*/
}

void MeshMaker::MakeParallelMesh() {
	pmesh = new ParMesh(MPI_COMM_WORLD, gmesh);
}

void MeshMaker::Make_H1_FESpace(int order) {
	H1_FECollection gFec(order, gmesh.Dimension());
	gFespace = new FiniteElementSpace(&gmesh, &gFec);
}	

void MeshMaker::Make_H1_FESpace_Parallel(int order) {
	H1_FECollection fec(order, pmesh->Dimension());
	fespace = new ParFiniteElementSpace(pmesh, &fec);
}
