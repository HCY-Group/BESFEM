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
	this->d->GetDerivative(1,0,gdX);
	this->d->GetDerivative(1,1,gdY);
		
	//calculate c ("velocity")
	ParGridFunction mgGd(this->d->ParFESpace());
	for (int vi = 0; vi < mgGd.Size(); vi++){
		// magnitude of gradient
		mgGd(vi) = sqrt( gdX(vi)*gdX(vi) + gdY(vi)*gdY(vi) );
		// c_x
		(*this->c)(vi) = (*this->sgn)(vi)*gdX(vi)/mgGd(vi);
		(*this->cx)(vi) = (*this->c)(vi);
		// c_y
		(*this->c)(vi+mgGd.Size()) = (*this->sgn)(vi)*gdY(vi)/mgGd(vi);
		(*this->cy)(vi) = (*this->c)(vi+mgGd.Size());
	}
		
}

void VoxelSolver_DG::FormMatrices() {

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

}
