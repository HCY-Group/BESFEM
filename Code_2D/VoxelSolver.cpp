
#include "VoxelSolver.hpp"
#include "TimeDepOpers.hpp"

using namespace mfem;
using namespace std;

VoxelSolver::VoxelSolver(FiniteElementSpace *fes){
	gVox = new GridFunction(fes);
	cout << "GVOX SIZE" << gVox->Size() << endl;
	
}
VoxelSolver::VoxelSolver(FiniteElementSpace *gfes, ParFiniteElementSpace *fes){
	gVox = new GridFunction(gfes);
	Vox = new ParGridFunction(fes);
	/*
	cout << "fespace: " << Vox->ParFESpace() << endl;
	cout << "fec:" << Vox->ParFESpace()->FEColl() << endl;
	cout << "finiteelement: " << Vox->ParFESpace()->GetFE(0) << endl;
	*/
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

void VoxelSolver::MapGlobalToLocal(Mesh* gmesh, ParMesh* pmesh) {
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
	
	GridFunction tmp_gf_glob(*gVox);	// Copy constructor (using to "dereference")
	ParGridFunction tmp_gf_par(*Vox);	// Copy constructor (using to "dereference")
	
	for (ei=0; ei<nE; ei++){
		gei = E_L2G[ei];

		gmesh->GetElementVertices(gei,gVTX);
		pmesh->GetElementVertices(ei,VTX);
	
		for (int vi = 0; vi<nC; vi++){
			//Vox(VTX[vi]) = gVox(gVTX[vi]);
			tmp_gf_par(VTX[vi]) = tmp_gf_glob(gVTX[vi]);
		}
	}
	
	*gVox = tmp_gf_glob;
	*Vox = tmp_gf_par;

}
/*
void VoxelSolver::InitStiffMatrix(Array<int> boundary_dofs, ParGridFunction Diff) {
	
	ParFiniteElementSpace fespace(*Diff.ParFESpace());
	
	ParLinearForm Fct(&fespace);
	//HypreParVector Fcb(&fespace);
	HypreParVector X1v(&fespace);
	
	Fcb = X1v.CreateCompatibleVector(); //needed so that Fcb is defined on a fespace?
	
	// stiffness matrix
	//HypreParMatrix Kmat;
	GridFunctionCoefficient cMob(&Diff);
	std::unique_ptr<ParBilinearForm> K(new ParBilinearForm(&fespace));
	K->AddDomainIntegrator(new DiffusionIntegrator(cMob));
	K->Assemble();
	K->FormLinearSystem(boundary_dofs, *Vox, Fct, Kmat, X1v, Fcb);
	
}

void VoxelSolver::InitTimeDepOper(ParGridFunction DomPar) {
	cout << "HERE A" << endl;
	// TimeDependentOperator and ODESolver
	cout << Fcb.Size() << endl;
	//ConductionOperator oper(DomPar, Kmat, Fcb);
	oper = new ConductionOperator(DomPar, Kmat, Fcb);
	//ODESolver *ode_solver = new ForwardEulerSolver;
	//ODESolver *ode_solver = new BackwardEulerSolver;
	ode_solver = new BackwardEulerSolver;
	ode_solver->Init(*oper);
	cout << "HERE B" << endl;
}
*/
void VoxelSolver::InitMatricesAndTimeDepOpers(Array<int> boundary_dofs, ParGridFunction Diff, ParGridFunction DomPar) {
	
	ParFiniteElementSpace fespace(*Diff.ParFESpace());
	cout << "fespace: " << Diff.ParFESpace() << endl;
	cout << "fec:" << Diff.ParFESpace()->FEColl() << endl;
	cout << "finiteelement: " << Diff.ParFESpace()->GetFE(0) << endl;
	
	//ParLinearForm Fct(&fespace);
	Fct = new ParLinearForm(&fespace);
	HypreParVector Fcb(&fespace);
	HypreParVector X1v(&fespace);
	
	Fcb = X1v.CreateCompatibleVector(); //needed so that Fcb is defined on a fespace?
	
	// stiffness matrix
	HypreParMatrix Kmat;
	GridFunctionCoefficient cMob(&Diff);
	//std::unique_ptr<ParBilinearForm> K(new ParBilinearForm(&fespace));
	K = new ParBilinearForm(&fespace);
	K->AddDomainIntegrator(new DiffusionIntegrator(cMob));
	K->Assemble();
	K->FormLinearSystem(boundary_dofs, *Vox, *Fct, Kmat, X1v, Fcb);

	//cout << "HERE A" << endl;
	// TimeDependentOperator and ODESolver
	cout << Fcb.Size() << endl;
	//ConductionOperator oper(DomPar, Kmat, Fcb);
	oper = new ConductionOperator(DomPar, Kmat, Fcb);
	//ODESolver *ode_solver = new ForwardEulerSolver;
	//ODESolver *ode_solver = new BackwardEulerSolver;
	ode_solver = new BackwardEulerSolver;
	ode_solver->Init(*oper);
	//cout << "HERE B" << endl;
}

void VoxelSolver::UpdateLinearForm(ParGridFunction gf) {
	
	GridFunctionCoefficient coef(&gf);
	//ParFiniteElementSpace fespace(*gf.ParFESpace());
	/*
	std::unique_ptr<ParLinearForm> Bc(new ParLinearForm(&fespace));
	Bc->AddDomainIntegrator(new DomainLFIntegrator(coef));
	Bc->Assemble();
	*/
	//ParLinearForm Bc(&fespace);
	ParLinearForm Bc(gf.ParFESpace());
	Bc.AddDomainIntegrator(new DomainLFIntegrator(coef));
	Bc.Assemble();
	//Fct = std::move(*Bc);
	Fct = std::move(&Bc);
}

void VoxelSolver::UpdateLinearForm_DoubleWellPotential() {
	///*
	cout << "HERE A" << endl;
	cout << "fespace: " << Vox->ParFESpace() << endl;
	cout << "fec:" << Vox->ParFESpace()->FEColl() << endl;
	cout << "finiteelement: " << Vox->ParFESpace()->GetFE(0) << endl;
	//cout << "finiteelement:" << Vox->ParFESpace()->FEColl()->GetFE() << endl;
	//ParFiniteElementSpace fespace(*Vox->ParFESpace());
	//ParFiniteElementSpace fespace(*fes_p);
	cout << "HERE B" << endl;
	//ParMesh *pmesh = fespace.GetParMesh();
	//int nV = pmesh->GetNV();
	//*/

	//ParGridFunction Pot(&fespace);
	ParGridFunction Pot(Vox->ParFESpace());
	ParGridFunction Vox2(*Vox);
	int nV = Pot.Size();
	
	for (int vi = 0; vi < nV; vi++){
		Pot(vi) = 2.0*Vox2(vi)*(1.0-Vox2(vi))*(1.0-2.0*Vox2(vi));
	}
	
	UpdateLinearForm(Pot);
	
}

void VoxelSolver::UpdateSystemAndSolve(Array<int> boundary_dofs, double t_ode, double dt) {
	
	ParFiniteElementSpace fespace(*Vox->ParFESpace());

	HypreParVector Fcb(&fespace);
	HypreParVector X1v(&fespace);
	HypreParVector Vox0(&fespace);
	HypreParMatrix Kmat;
	
	Vox->GetTrueDofs(Vox0);
	K->FormLinearSystem(boundary_dofs, *Vox, *Fct, Kmat, X1v, Fcb);
	oper->UpdateParams(Kmat, Fcb);
	ode_solver->Step(Vox0, t_ode, dt);
	Vox->Distribute(Vox0);
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
