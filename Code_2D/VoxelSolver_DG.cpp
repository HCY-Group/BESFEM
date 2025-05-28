#include "VoxelSolver_DG.hpp"


VoxelSolver_DG::VoxelSolver_DG(FiniteElementSpace *gfes, ParFiniteElementSpace *fes)
	: VoxelSolver(gfes, fes){

}

VoxelSolver_DG::VoxelSolver_DG(FiniteElementSpace *gfes, ParFiniteElementSpace *fes_h1, ParFiniteElementSpace *fes_dg, ParFiniteElementSpace *dimfes_dg)
	: VoxelSolver(gfes, fes_h1){

	this->d   = new ParGridFunction(fes_dg);
	this->c   = new ParGridFunction(dimfes_dg);
	this->cx  = new ParGridFunction(fes_dg);
	this->cy  = new ParGridFunction(fes_dg);
	//cout << "VECTORDIM: " << this->c->VectorDim() << endl;
	if ( this->c->VectorDim() == 3 ){
		this->cz  = new ParGridFunction(fes_dg);
	}
	this->sgn = new ParGridFunction(fes_dg);
}

void VoxelSolver_DG::ProjectVals(ParGridFunction *gf) {
	
	this->d->ProjectGridFunction(*gf);
	
	for (int vi = 0; vi < this->sgn->Size(); vi++){
		if ((*this->d)(vi)<0.5){
			(*this->sgn)(vi) = -1;
		} else {
			(*this->sgn)(vi) = 1;
		}
	}
}

void VoxelSolver_DG::CalcLevelSetVel() {
	
	//gradient of d
	ParGridFunction gdX(this->d->ParFESpace());
	ParGridFunction gdY(this->d->ParFESpace());
	ParGridFunction gdZ(this->d->ParFESpace());
	this->d->GetDerivative(1,0,gdX);
	this->d->GetDerivative(1,1,gdY);
	if ( this->c->VectorDim() == 3 ){
		this->d->GetDerivative(1,2,gdZ);
	}
	//calculate c ("velocity")
	ParGridFunction mgGd(this->d->ParFESpace());
	for (int vi = 0; vi < mgGd.Size(); vi++){
		// magnitude of gradient
		if ( this->c->VectorDim() ==  2 ){
			mgGd(vi) = sqrt( gdX(vi)*gdX(vi) + gdY(vi)*gdY(vi) );
		}
		if ( this->c->VectorDim() ==  3 ){
			mgGd(vi) = sqrt( gdX(vi)*gdX(vi) + gdY(vi)*gdY(vi) + gdZ(vi)*gdZ(vi) );
		}
		// c_x
		(*this->cx)(vi) = (*this->sgn)(vi)*gdX(vi)/mgGd(vi);
		// c_y
		(*this->cy)(vi) = (*this->sgn)(vi)*gdY(vi)/mgGd(vi);
		if ( this->c->VectorDim() ==  3 ){
			(*this->cz)(vi) = (*this->sgn)(vi)*gdZ(vi)/mgGd(vi);
		}
		if (mgGd(vi) < 1e-3){
			(*this->cx)(vi) = 0.0;
			(*this->cy)(vi) = 0.0;
			if ( this->c->VectorDim() ==  3 ){
				(*this->cz)(vi) = 0.0;
			}
		}
		(*this->c)(vi) = (*this->cx)(vi);
		(*this->c)(vi+mgGd.Size()) = (*this->cy)(vi);
		if ( this->c->VectorDim() ==  3 ){
			(*this->c)(vi+2*mgGd.Size()) = (*this->cz)(vi);
		}
	}
		
}

void VoxelSolver_DG::FormMatrices(Array<int> ess_tdof_list) {

	//calculate M matrix
	this->m = new ParBilinearForm(this->d->ParFESpace());
	this->m->AddDomainIntegrator(new MassIntegrator);
	
	//calculate K matrix
	real_t alpha = -1.0;
	this->k = new ParBilinearForm(this->d->ParFESpace());
	VectorGridFunctionCoefficient cCoef(this->c);
   	this->k->AddDomainIntegrator(new ConvectionIntegrator(cCoef, alpha));
  	this->k->AddInteriorFaceIntegrator(
      		new NonconservativeDGTraceIntegrator(cCoef, alpha));
		
	//Add a (small) diffusive term to allow Neumann BCs to be enforced weakly by FEM
	ConstantCoefficient diffCoef(1e-8);
	this->k->AddDomainIntegrator(new DiffusionIntegrator(diffCoef));
	this->k->AddBdrFaceIntegrator(
		new DGDiffusionIntegrator(diffCoef,-1,-1));
	this->k->AddInteriorFaceIntegrator(
		new DGDiffusionIntegrator(diffCoef,-1,-1));

	//Form b vector
	this->b = new ParLinearForm(this->d->ParFESpace());
	GridFunctionCoefficient sgnCoef(this->sgn);
	this->b->AddDomainIntegrator(new DomainLFIntegrator(sgnCoef));

	//Assemble matrices and vectors
   	int skip_zeros = 0;
   	this->m->Assemble();
   	this->k->Assemble(skip_zeros);
   	this->b->Assemble();
   	this->m->Finalize();
   	this->k->Finalize(skip_zeros);

   	HypreParVector *B = this->b->ParallelAssemble();
	HypreParMatrix Kmat;
	this->k->FormSystemMatrix(ess_tdof_list, Kmat);
	this->advec = new AdvectionOperator(*this->d, Kmat, *B);
	
	this->ode_solver = new ForwardEulerSolver;
	ode_solver->Init(*advec);
}


void VoxelSolver_DG::UpdateMatricesAndSolve(Array<int> ess_tdof_list, double t_ode, double dt) {

	//calculate K matrix
	real_t alpha = -1.0;
	delete this->k;
	this->k = new ParBilinearForm(this->d->ParFESpace());
	VectorGridFunctionCoefficient cCoef(this->c);
   	this->k->AddDomainIntegrator(new ConvectionIntegrator(cCoef, alpha));
  	this->k->AddInteriorFaceIntegrator(
      		new NonconservativeDGTraceIntegrator(cCoef, alpha));
		
	//Add a (small) diffusive term to allow Neumann BCs to be enforced weakly by FEM
	ConstantCoefficient diffCoef(1e-8);
	this->k->AddDomainIntegrator(new DiffusionIntegrator(diffCoef));
	this->k->AddBdrFaceIntegrator(
		new DGDiffusionIntegrator(diffCoef,-1,-1));
	this->k->AddInteriorFaceIntegrator(
		new DGDiffusionIntegrator(diffCoef,-1,-1));

	//Form b vector
	delete this->b;
	this->b = new ParLinearForm(this->d->ParFESpace());
	GridFunctionCoefficient sgnCoef(this->sgn);
	this->b->AddDomainIntegrator(new DomainLFIntegrator(sgnCoef));

	//Assemble matrices and vectors
   	int skip_zeros = 0;
   	this->k->Assemble(skip_zeros);
   	this->b->Assemble();
   	this->k->Finalize(skip_zeros);
	
	/*
   	HypreParVector *B = this->b->ParallelAssemble();
	HypreParMatrix Kmat;
	this->k->FormSystemMatrix(ess_tdof_list, Kmat);
	
	this->advec->UpdateParams(Kmat,*B);
	*/
	HypreParVector Fcb(this->d->ParFESpace());
	HypreParVector X1v(this->d->ParFESpace());
	HypreParMatrix Kmat;
	
	this->k->FormLinearSystem(ess_tdof_list, *this->d, *this->b, Kmat, X1v, Fcb);
	this->d->GetTrueDofs(X1v); //Set initial values
	advec->UpdateParams(Kmat, Fcb);
	ode_solver->Step(X1v, t_ode, dt);
	//cout << "Max d:  " << X1v.Max() << ", Min d:  " << X1v.Min() << endl;
	this->k->RecoverFEMSolution(X1v, *this->b, *this->d);
	//ode_solver->Step(d, t_ode, dt);
	//cout << "Max d:  " << d->Max() << ", Min d:  " << d->Min() << endl;

}
