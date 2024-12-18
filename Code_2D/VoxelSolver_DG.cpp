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
