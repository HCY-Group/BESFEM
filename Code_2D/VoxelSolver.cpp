
#include "VoxelSolver.hpp"
#include "TimeDepOpers.hpp"

using namespace mfem;
using namespace std;

VoxelSolver::VoxelSolver(mfem::FiniteElementSpace *fes){
	this->gVox = new mfem::GridFunction(fes);
	
}
VoxelSolver::VoxelSolver(mfem::FiniteElementSpace *gfes, mfem::ParFiniteElementSpace *fes){
	this->gVox = new mfem::GridFunction(gfes);
	this->Vox = new mfem::ParGridFunction(fes);
}
	
void VoxelSolver::AssignGlobalValues(double value){
	*this->gVox = value;
}

void VoxelSolver::AssignGlobalValues(vector<vector<vector<int>>> data) {

	// global grid function for voxel data
	cout << "Defining Voxel GridFunction" << endl;
	int nz = data.size();
	int ny = data[0].size();
	int nx = data[0][0].size();
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				int idx = i + nx*j + nx*ny*k;
				(*this->gVox)[idx] = data[k][j][i];
			}
		}
	}
}

void VoxelSolver::MultiplyVox(mfem::GridFunction &gf) { 
	//Elementwise Multiplication
	*this->Vox *= gf;
}

void VoxelSolver::MapGlobalToLocal(mfem::Mesh* gmesh, mfem::ParMesh* pmesh) {
	int nV = pmesh->GetNV();		//number of vertices
	int nE = pmesh->GetNE();		//number of elements
	int nC = pow(2,pmesh->Dimension());	//number of corner vertices
	
	// Map local to global element indices
	Array<HYPRE_BigInt> E_L2G;
	pmesh->GetGlobalElementIndices(E_L2G);

	Array<int> gVTX(nC);	//global indices of corner vertices
	Array<int> VTX(nC);	//local indices of corner vertices
	int gei;		//global element indices
	int ei;			//local element indices
	
	
	for (ei=0; ei<nE; ei++){
		gei = E_L2G[ei];

		gmesh->GetElementVertices(gei,gVTX);
		pmesh->GetElementVertices(ei,VTX);
	
		for (int vi = 0; vi<nC; vi++){
			(*this->Vox)(VTX[vi]) = (*this->gVox)(gVTX[vi]);
		}
	}
	

}

void VoxelSolver::AssignDirichletBCs(mfem::GridFunctionCoefficient &Coef, Array<int> &bdr) {
	this->Vox->ProjectBdrCoefficient(Coef, bdr);
}

void VoxelSolver::DetermineConnectivityBCs(mfem::ParGridFunction &DomPar) {
	if (this->dbcval != nullptr) {delete this->dbcval;}
	this->dbcval = new mfem::ParGridFunction(this->Vox->ParFESpace());
	*this->dbcval = 0.0;

	for (int vi = 0; vi < this->Vox->Size(); vi++) {
		if ( DomPar(vi) > 0.6 ) {
			(*this->dbcval)(vi) = 1.0;
		}
	}
	
	if (this->dbcCoef != nullptr) {delete this->dbcCoef;}
	this->dbcCoef = new mfem::GridFunctionCoefficient(this->dbcval);
}

void VoxelSolver::NorthDirichletBCs(mfem::Mesh *mesh) {
	if (this->dbc_bdr != nullptr) {delete this->dbc_bdr;}
	this->dbc_bdr = new Array<int>;
	this->dbc_bdr->SetSize(mesh->bdr_attributes.Max());
	*this->dbc_bdr = 0;
	cout << "HERE A" << endl;
	cout << "Size: " << this->dbc_bdr->Size() << endl;
	if (mesh->Dimension()==2) {
		(*this->dbc_bdr)[2] = 1;
	} else {
		(*this->dbc_bdr)[3] = 1;
	}
	cout << "HERE B" << endl;

	if (this->ess_tdof_list != nullptr) {delete this->ess_tdof_list;}
	this->ess_tdof_list = new Array<int>;
	cout << "HERE C" << endl;
	this->Vox->ParFESpace()->GetEssentialTrueDofs(*this->dbc_bdr, *this->ess_tdof_list);
	
}

void VoxelSolver::SouthDirichletBCs(mfem::Mesh *mesh) {
	if (this->dbc_bdr != nullptr) {delete this->dbc_bdr;}
	this->dbc_bdr = new Array<int>;
	this->dbc_bdr->SetSize(mesh->bdr_attributes.Max());
	*this->dbc_bdr = 0;
	cout << "HERE A" << endl;
	cout << "Size: " << this->dbc_bdr->Size() << endl;
	if (mesh->Dimension()==2) {
		(*this->dbc_bdr)[0] = 1;
	} else {
		(*this->dbc_bdr)[1] = 1;
	}
	cout << "HERE B" << endl;

	if (this->ess_tdof_list != nullptr) {delete this->ess_tdof_list;}
	this->ess_tdof_list = new Array<int>;
	cout << "HERE C" << endl;
	this->Vox->ParFESpace()->GetEssentialTrueDofs(*this->dbc_bdr, *this->ess_tdof_list);
	
}

void VoxelSolver::InitMatricesAndTimeDepOpers(Array<int> boundary_dofs, mfem::ParGridFunction &Diff, mfem::ParGridFunction &DomPar) {

	// TODO: Add check to make sure that Diff and DomPar have same FESpace as Vox
	this->Fct = new mfem::ParLinearForm(Vox->ParFESpace());
	mfem::HypreParVector Fcb(Vox->ParFESpace());
	mfem::HypreParVector X1v(Vox->ParFESpace());
	

	Fcb = X1v.CreateCompatibleVector(); //needed so that Fcb is defined on a fespace?
	
	// stiffness matrix
	mfem::HypreParMatrix Kmat;
	mfem::GridFunctionCoefficient cMob(&Diff);
	K = new mfem::ParBilinearForm(Diff.ParFESpace());
	K->AddDomainIntegrator(new mfem::DiffusionIntegrator(cMob));
	K->Assemble();
	K->FormLinearSystem(boundary_dofs, *Vox, *Fct, Kmat, X1v, Fcb);

	// TimeDependentOperator and ODESolver
	cout << Fcb.Size() << endl;
	oper = new ConductionOperator(DomPar, Kmat, Fcb, boundary_dofs);
	ode_solver = new mfem::BackwardEulerSolver;
	ode_solver->Init(*oper);
}

void VoxelSolver::UpdateLinearForm(mfem::ParGridFunction gf) {
	
	mfem::GridFunctionCoefficient coef(&gf);
	
	delete this->Fct;
	this->Fct = new mfem::ParLinearForm(this->Vox->ParFESpace());
	this->Fct->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));
	this->Fct->Assemble();
	
}

void VoxelSolver::UpdateLinearForm_DoubleWellPotential() {

	mfem::ParGridFunction Pot(Vox->ParFESpace());
	int nV = Pot.Size();
	
	for (int vi = 0; vi < nV; vi++){
		Pot(vi) = 2.0*(*Vox)(vi)*(1.0-(*Vox)(vi))*(1.0-2.0*(*Vox)(vi));
	}
	
	UpdateLinearForm(Pot);
	
}

void VoxelSolver::UpdateSystemAndSolve(Array<int> boundary_dofs, double t_ode, double dt) {
	
	mfem::HypreParVector Fcb(this->Vox->ParFESpace());
	mfem::HypreParVector X1v(this->Vox->ParFESpace());
	mfem::HypreParVector Vox0(this->Vox->ParFESpace());
	mfem::HypreParMatrix Kmat;
	
	//cout << "Max Fct: " << Fct->Max() << ", Min Fct:  " << Fct->Min() << endl;
	K->FormLinearSystem(boundary_dofs, *Vox, *Fct, Kmat, X1v, Fcb);
	Vox->GetTrueDofs(X1v); //Set initial values
	//cout << "Size Fct: " << Fct->Size() << ", Size Vox0: " << Vox0.Size() << endl;
	oper->UpdateParams(Kmat, Fcb);
	ode_solver->Step(X1v, t_ode, dt);
	cout << "Max Vox0: " << X1v.Max() << ", Min Vox0:  " << X1v.Min() << endl;
	K->RecoverFEMSolution(X1v, *Fct, *Vox);
}

void VoxelSolver::AccelerateDiffusion(mfem::ParGridFunction &DomPar, mfem::GridFunctionCoefficient &Coef, Array<int> &bdr) {
	//TODO: add check to see if DomPar and Vox are same size
	
	int nV = Vox->Size();
	for (int vi = 0; vi < nV; vi++){
		if ( (*Vox)(vi) > 0.1 && DomPar(vi) > 0.9 ){
			(*Vox)(vi) = 1.0; // Modify Vox
		}
	}
	
	AssignDirichletBCs(Coef, bdr);
}

void VoxelSolver::ParaviewSave(string FileName, string VariableName, mfem::GridFunction* gf) {
	
	int order = 1;
	
	mfem::ParaViewDataCollection *pd = NULL;
	pd = new mfem::ParaViewDataCollection(FileName, gf->FESpace()->GetMesh());
	pd->RegisterField(VariableName, gf);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;
}
