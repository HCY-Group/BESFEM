
#include "MeshMaker.hpp"

using namespace mfem;
using namespace std;

MeshMaker::MeshMaker(std::vector<std::vector<std::vector<int>>> tiffdata){
	this->data = tiffdata;
}

void MeshMaker::MakeGlobalMesh() {
	int nz = this->data.size();
	int ny = this->data[0].size();
	int nx = this->data[0][0].size();
	double sx = nx;  //make dx = 1
	double sy = ny;  //make dy = 1
	double sz = nz;  //make dz = 1
	bool generate_edges = false;
	bool sfc_ordering = false;
	this->gmesh = new Mesh;
	if (nz == 1) {
		*this->gmesh = Mesh::MakeCartesian2D(nx-1, ny-1, Element::QUADRILATERAL, generate_edges, sx, sy, sfc_ordering);
	} else {
		*this->gmesh = Mesh::MakeCartesian3D(nx-1, ny-1, nz-1, Element::HEXAHEDRON, sx, sy, sz, sfc_ordering);
	}
	this->gmesh->EnsureNCMesh(true);

	cout << "Mesh from MeshMaker " << this->gmesh << endl;
}

void MeshMaker::MakeParallelMesh() {
	this->pmesh = new ParMesh(MPI_COMM_WORLD, *this->gmesh);
}

void MeshMaker::Make_H1_FESpace(int order) {
	H1_FECollection gFec(order, this->gmesh->Dimension());
	this->gFespace = new FiniteElementSpace(this->gmesh, &gFec);
}	

void MeshMaker::Make_H1_FESpace_Parallel(int order) {
	H1_FECollection fec(order, this->pmesh->Dimension());
	this->fespace = new ParFiniteElementSpace(this->pmesh, &fec);
	
	cout << "fespace: " << this->fespace << endl;
	cout << "pmesh: " << this->fespace->GetParMesh() << endl;
	cout << "fec: " << this->fespace->FEColl() << endl;
	
	cout << "finiteelement: " << this->fespace->GetFE(0) << endl;
}
