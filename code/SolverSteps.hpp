#ifndef SOLVERSTEPS_HPP
#define SOLVERSTEPS_HPP

#include "Initialize_Geometry.hpp"
#include "Constants.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory>

class SolverSteps {
public:
    SolverSteps(std::shared_ptr<mfem::ParFiniteElementSpace> fespace);
    ~SolverSteps();

    void MassMatrix(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat);
    
    void StiffnessMatrix(std::shared_ptr<mfem::GridFunctionCoefficient> cDx, mfem::Array<int> boundary, 
                         mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, 
                         std::shared_ptr<mfem::HypreParMatrix> &Kmat, mfem::HypreParVector &RHS);
    
    void ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F);
    void ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, mfem::Array<int> &boundary, mfem::ProductCoefficient &m);
    
    void SolverConditions(std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &solver, 
                          mfem::HypreSmoother &smoother);

    std::shared_ptr<mfem::ParBilinearForm> K;

private:
    mfem::ParMesh *pmesh;
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;
        
    mfem::ParBilinearForm *M = nullptr;
    std::unique_ptr<mfem::ParLinearForm> B;
    mfem::ParGridFunction *temp_ps = nullptr;
    mfem::GridFunctionCoefficient *coef = nullptr;
    mfem::HypreParVector X1v;
    mfem::Array<int> boundary_dofs; ///< Boundary degrees of freedom

};

#endif // SOLVER_HPP
