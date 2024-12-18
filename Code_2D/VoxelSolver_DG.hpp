#ifndef VOXELSOLVER_DG_HPP
#define VOXELSOLVER_DG_HPP

#include "VoxelSolver.hpp"

class VoxelSolver_DG : public VoxelSolver {
public:
	VoxelSolver_DG(FiniteElementSpace *gfes, ParFiniteElementSpace *fes); //constructor with parallel
	VoxelSolver_DG(FiniteElementSpace *gfes, ParFiniteElementSpace *fes_h1, ParFiniteElementSpace *fes_dg, ParFiniteElementSpace *dimfes_dg); //constructor with parallel
	void ProjectVals(ParGridFunction *gf);
	void CalcLevelSetVel();
	void FormMatrices();

private:
	ParGridFunction *d = nullptr;
	ParGridFunction *c = nullptr;
	ParGridFunction *cx = nullptr;
	ParGridFunction *cy = nullptr;
	ParGridFunction *sgn = nullptr;

	ParBilinearForm *m = nullptr;
	ParBilinearForm *k = nullptr;
	ParLinearForm *b = nullptr;
};

#endif
