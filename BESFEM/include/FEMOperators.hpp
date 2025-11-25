#ifndef FEMOPERATORS_HPP
#define FEMOPERATORS_HPP

#include "Initialize_Geometry.hpp"
#include "Constants.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory>

/**
 * @class FEMOperators
 * @brief Utility class for constructing and manipulating FEM operators in BESFEM.
 *
 * Provides reusable routines for:
 * - Mass and stiffness matrix assembly (ParBilinearForm)
 * - Force/forcing term assembly (ParLinearForm)
 * - Forming system matrices with essential boundary conditions
 * - Assembling linear systems (A, X, B) for Hypre solvers
 * - Setting solver and preconditioner options
 *
 * This class centralizes common MFEM assembly tasks used by CnA, CnC, and CnE
 * to reduce duplication across the concentration solvers.
 */
class FEMOperators {
public:

    /**
     * @brief Construct FEMOperators using a shared FE space.
     *
     * @param fespace Parallel finite element space (shared_ptr).
     */
    FEMOperators(std::shared_ptr<mfem::ParFiniteElementSpace> fespace);

    /**
     * @brief Construct FEMOperators from a raw FE space pointer.
     *
     * @param fespace Raw pointer to a parallel FE space.
     */
    FEMOperators(mfem::ParFiniteElementSpace *fespace);

    /// Destructor.
    ~FEMOperators();

    /**
     * @brief Initialize a mass bilinear form.
     *
     * Creates a ParBilinearForm on the FE space using coefficient @p coef and
     * performs domain integration to build the mass matrix structure.
     *
     * @param coef Scalar or spatial coefficient.
     * @param M    Output unique_ptr to the mass bilinear form.
     */
    void InitializeMassMatrix(mfem::Coefficient &coef, std::unique_ptr<mfem::ParBilinearForm> &M);

    /**
     * @brief Initialize a stiffness bilinear form.
     *
     * Creates and assembles a diffusion/gradient operator using the provided
     * coefficient @p coef.
     *
     * @param coef Coefficient used for the stiffness matrix.
     * @param K    Output unique_ptr to the stiffness bilinear form.
     */
    void InitializeStiffnessMatrix(mfem::Coefficient &coef, std::unique_ptr<mfem::ParBilinearForm> &K);

    /**
     * @brief Initialize a force (linear) form.
     *
     * Constructs and assembles a ParLinearForm using coefficient @p coef.
     *
     * @param coef Coefficient for the force term.
     * @param B    Output unique_ptr to the linear form.
     */
    void InitializeForceTerm(mfem::Coefficient &coef, std::unique_ptr<mfem::ParLinearForm> &B);

    /**
     * @brief Form the Hypre system matrix from a bilinear form.
     *
     * Applies essential boundary DOFs, enforcing homogeneous Dirichlet
     * conditions, and exports the final system matrix into @p A.
     *
     * @param M             Mass/stiffness bilinear form.
     * @param boundary_dofs List of essential DOFs.
     * @param A             Output HypreParMatrix system matrix.
     */
    void FormSystemMatrix(std::unique_ptr<mfem::ParBilinearForm> &M, mfem::Array<int> &boundary_dofs, mfem::HypreParMatrix &A);

    /**
     * @brief Form the linear system A·X = B.
     *
     * Combines stiffness/mass matrices with boundary conditions, extracts the
     * Hypre system matrix and RHS vector, and maps the solution grid-function
     * @p x into true DOFs.
     *
     * @param K             Bilinear form (stiffness or mass).
     * @param boundary_dofs Essential DOFs.
     * @param x             Solution grid-function.
     * @param b             Linear form (load/forcing).
     * @param A             Output system matrix.
     * @param X             Output solution vector (true DOFs).
     * @param B             Output RHS vector (true DOFs).
     */
    void FormLinearSystem(std::unique_ptr<mfem::ParBilinearForm> &K, mfem::Array<int> &boundary_dofs, mfem::ParGridFunction &x,
                          mfem::ParLinearForm &b, mfem::HypreParMatrix &A, mfem::HypreParVector &X, mfem::HypreParVector &B);

    /**
     * @brief Configure solver and preconditioner for a given system matrix.
     *
     * Sets tolerances, iteration caps, and attaches @p preconditioner to @p solver.
     *
     * @param Mmat          System matrix.
     * @param solver        Hypre CG solver to configure.
     * @param preconditioner Preconditioner (Jacobi, Gauss-Seidel, HypreSmoother, etc.).
     */
    void SolverConditions(mfem::HypreParMatrix &Mmat, mfem::CGSolver &solver, mfem::Solver &preconditioner);

    /**
     * @brief Configure solver/preconditioner pair (matrix-free).
     *
     * Useful when matrix is already installed in the solver.
     *
     * @param solver        CG solver.
     * @param preconditioner Preconditioner.
     */
    void SolverConditions(mfem::CGSolver &solver, mfem::Solver &preconditioner);

    /**
     * @brief Reassemble a ParLinearForm.
     *
     * @param B Linear form to update.
     */
    void Update(std::unique_ptr<mfem::ParLinearForm> &B);

    /**
     * @brief Reassemble a ParBilinearForm.
     *
     * @param B Bilinear form to update.
     */
    void Update(std::unique_ptr<mfem::ParBilinearForm> &B);

    // -------------------------------------------------------------------------
    // Public fields
    // -------------------------------------------------------------------------
    std::shared_ptr<mfem::ParBilinearForm> K; ///< Assembled stiffness form.
    mfem::HypreParVector X1v;                 ///< Scratch true-DoF vector.
    mfem::Array<int> boundary_dofs;           ///< Essential boundary DOFs.

private:
    // -------------------------------------------------------------------------
    // Internal FE/Mesh context
    // -------------------------------------------------------------------------
    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< FE space (shared_ptr).
    mfem::ParFiniteElementSpace *raw_fespace = nullptr;   ///< Raw FE space pointer.
    mfem::ParFiniteElementSpace *local_fespace = nullptr; ///< Local FE space (alias).

    // -------------------------------------------------------------------------
    // Operator storage
    // -------------------------------------------------------------------------
    mfem::ParBilinearForm *M = nullptr; ///< Mass matrix (raw pointer).
    std::unique_ptr<mfem::ParLinearForm> B; ///< Force/forcing term.
    mfem::ParGridFunction *temp_ps = nullptr; ///< Temporary phase-field grid-function.
    mfem::GridFunctionCoefficient *coef = nullptr; ///< Coefficient wrapper for grid-functions.
};

#endif // FEMOPERATORS_HPP
