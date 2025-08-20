#ifndef SOLVERSTEPS_HPP
#define SOLVERSTEPS_HPP

/**
 * @file SolverSteps.hpp
 * @brief Base class for common solver steps in battery simulations.
 *
 * Provides methods for initializing mass/stiffness matrices, assembling
 * linear systems, and applying solvers/preconditioners.
 */

#include "Initialize_Geometry.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory>

/**
 * @class SolverSteps
 * @brief Base class for common solver steps in battery simulations.
 *
 * Provides methods for initializing mass/stiffness matrices, assembling
 * linear systems, and applying solvers/preconditioners.
 */
class SolverSteps {
public:

    /**
     * @brief Constructor for SolverSteps
     * @param fespace Pointer to the finite element space
     */
    SolverSteps(std::shared_ptr<mfem::ParFiniteElementSpace> fespace);
    
    /**
     * @brief Constructor for SolverSteps
     * @param fespace Pointer to the finite element space
     */
    SolverSteps(mfem::ParFiniteElementSpace* fespace);
    
    ~SolverSteps();

    /**
     * @brief Initializes the mass matrix.
     * @param coef Coefficient for the mass matrix.
     * @param M Pointer to the mass matrix.
     */
    void InitializeMassMatrix(mfem::Coefficient &coef, std::unique_ptr<mfem::ParBilinearForm> &M);
    

    /**
     * @brief Initializes the stiffness matrix.
     * @param coef Coefficient for the stiffness matrix.
     * @param K Pointer to the stiffness matrix.
     */
    void InitializeStiffnessMatrix(mfem::Coefficient &coef, std::unique_ptr<mfem::ParBilinearForm> &K);
    
    
    /**
     * @brief Initializes the force term.
     * @param coef Coefficient for the force term.
     * @param B Pointer to the linear form for the force term.
     * @param boundary_coef Optional boundary coefficient.
     * @param boundary_attr Optional boundary attributes.
     */
    void InitializeForceTerm(mfem::Coefficient &coef, std::unique_ptr<mfem::ParLinearForm> &B, mfem::Coefficient *boundary_coef = nullptr,
        mfem::Array<int> *boundary_attr = nullptr);
    
    /**
     * @brief Forms the system matrix.
     * @param M Pointer to the mass matrix.
     * @param boundary_dofs Boundary degrees of freedom.
     * @param A Output system matrix.
     */
    void FormSystemMatrix(std::unique_ptr<mfem::ParBilinearForm> &M, mfem::Array<int> &boundary_dofs, mfem::HypreParMatrix &A);
    
    /**
     * @brief Forms the linear system.
     * @param K Pointer to the stiffness matrix.
     * @param boundary_dofs Boundary degrees of freedom.
     * @param x Solution grid function.
     * @param b Right-hand side linear form.
     * @param A Output system matrix.
     * @param X Output solution vector.
     * @param B Output right-hand side vector.
     */
    void FormLinearSystem(std::unique_ptr<mfem::ParBilinearForm> &K, mfem::Array<int> &boundary_dofs, mfem::ParGridFunction &x, mfem::ParLinearForm &b, mfem::HypreParMatrix &A, mfem::HypreParVector &X, mfem::HypreParVector &B);
    
    /**
     * @brief Sets up solver conditions.
     * @param Mmat System matrix.
     * @param solver Solver to configure.
     * @param preconditioner Preconditioner to configure.
     */
    void SolverConditions(mfem::HypreParMatrix &Mmat, mfem::CGSolver &solver, mfem::Solver &preconditioner);
    
    /**
     * @brief Sets up solver conditions.
     * @param solver Solver to configure.
     * @param preconditioner Preconditioner to configure.
     */
    void SolverConditions(mfem::CGSolver &solver, mfem::Solver &preconditioner);

    /**
     * @brief Updates the Parlinear form.
     * @param B Parlinear form to update.
     */
    void Update(std::unique_ptr<mfem::ParLinearForm> &B);

    /**
     * @brief Updates the ParBilinear form.
     * @param B ParBilinear form to update.
     */
    void Update(std::unique_ptr<mfem::ParBilinearForm> &B);
    

    std::shared_ptr<mfem::ParBilinearForm> K; ///< Pointer to the stiffness matrix
    mfem::HypreParVector X1v; ///< Solution vector

private:
    mfem::ParMesh *pmesh; ///< Pointer to the parallel mesh
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space
    mfem::ParFiniteElementSpace* raw_fespace; ///< Pointer to the raw finite element space
    mfem::ParFiniteElementSpace* local_fespace; ///< Pointer to the local finite element space
        
    mfem::ParBilinearForm *M = nullptr; ///< Pointer to the mass matrix
    std::unique_ptr<mfem::ParLinearForm> B; ///< Pointer to the linear form for the force term
    mfem::ParGridFunction *temp_ps = nullptr; ///< Pointer to a temporary grid function for phase field
    mfem::GridFunctionCoefficient *coef = nullptr; ///< Pointer to a coefficient for the grid function
    mfem::Array<int> boundary_dofs; ///< Boundary degrees of freedom
};

#endif // SOLVERSTEPS_HPP