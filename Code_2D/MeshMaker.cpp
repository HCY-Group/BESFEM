
#include "MeshMaker.hpp"


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
	this->gmesh = new mfem::Mesh;
	if (nz == 1) {
		*this->gmesh = mfem::Mesh::MakeCartesian2D(nx-1, ny-1, mfem::Element::QUADRILATERAL, generate_edges, sx, sy, sfc_ordering);
	} else {
		*this->gmesh = mfem::Mesh::MakeCartesian3D(nx-1, ny-1, nz-1, mfem::Element::HEXAHEDRON, sx, sy, sz, sfc_ordering);
	}
	this->gmesh->EnsureNCMesh(true);

}

void MeshMaker::MakeParallelMesh() {
	this->pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *this->gmesh);
}

void MeshMaker::Make_H1_FESpace(int order) {
	this->gFec = new mfem::H1_FECollection(order, this->gmesh->Dimension());
	this->gFespace = new mfem::FiniteElementSpace(this->gmesh, this->gFec);
}	

void MeshMaker::Make_H1_FESpace_Parallel(int order) {
	this->fec = new mfem::H1_FECollection(order, this->pmesh->Dimension());
	this->fespace = new mfem::ParFiniteElementSpace(this->pmesh, this->fec);
}

void MeshMaker::Make_DG_FESpace_Parallel(int order) {
	
	this->fec_dg = new mfem::DG_FECollection(order, this->pmesh->Dimension(), mfem::BasisType::GaussLobatto);
	this->fespace_dg = new mfem::ParFiniteElementSpace(this->pmesh, this->fec_dg);
	this->dimfespace_dg = new mfem::ParFiniteElementSpace(this->pmesh, this->fec_dg, this->pmesh->Dimension(), mfem::Ordering::byNODES);
	
}

void MeshMaker::TestGetFE() {
	std::cout << "TestGetFE" << std::endl;
	std::cout << "fespace: " << this->fespace << std::endl;
	std::cout << "pmesh: " << this->fespace->GetParMesh() << std::endl;
	std::cout << "fec: " << this->fespace->FEColl() << std::endl;
	std::cout << "finiteelement: " << this->fespace->GetFE(0) << std::endl;
}







