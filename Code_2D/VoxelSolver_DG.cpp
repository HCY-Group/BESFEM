#include "VoxelSolver_DG.hpp"


VoxelSolver_DG::VoxelSolver_DG(FiniteElementSpace *gfes, ParFiniteElementSpace *fes)
	: VoxelSolver(gfes, fes){

}

VoxelSolver_DG::VoxelSolver_DG(FiniteElementSpace *gfes, ParFiniteElementSpace *fes_h1, ParFiniteElementSpace *fes_dg, ParFiniteElementSpace *dimfes_dg)
	: VoxelSolver(gfes, fes_h1){

	this->d  = new ParGridFunction(fes_dg);
	this->c  = new ParGridFunction(dimfes_dg);
	this->cx = new ParGridFunction(fes_dg);
	this->cy = new ParGridFunction(fes_dg);
}
