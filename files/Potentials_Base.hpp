#ifndef POTENTIALS_HPP
#define POTENTIALS_HPP

// Public: Members can be accessed from anywhere. This is the default access modifier. 
// Protected: Members can be accessed within the class and by classes that inherit from that class. 
// Private: Members can only be accessed within the class that defines them.


#include "mfem.hpp"
#include "Mesh_Handler.hpp"

#include <memory>

/**
 * @class Potentials
 * @brief Base class for handling potentials in battery simulations.
 *
 * This class defines methods to initialize potentials, set up solvers,
 * and perform operations like reaction generation, force term assembly,
 * and error calculations for battery simulation.
 */
class Potentials {
public:

    /**
     * @brief Constructor for the Potentials class.
     * @param pmesh Pointer to the parallel mesh.
     * @param fespace Pointer to the finite element space.
     * @param mh Reference to the MeshHandler for geometry and volume data.
     */
    Potentials(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);

    virtual ~Potentials() = default; ///< Virtual destructor.

    MeshHandler &mesh_handler; ///< Reference to the MeshHandler for geometry and volume operations.

    /**
     * @brief Initializes the potentials to a specified initial value.
     * @param ph Grid function representing the potential field.
     * @param initial_value Initial value for the potentials.
     */
    void SetInitialPotentials(mfem::ParGridFunction &ph, double initial_value);
    
    /**
     * @brief Configures the solver with given tolerance and iteration limits.
     * @param solver Conjugate gradient solver instance.
     * @param value_1 Relative tolerance for convergence.
     * @param value_2 Maximum number of iterations.
     */
    void SetUpSolver(mfem::CGSolver &solver, double value_1, double value_2);
    
    /**
     * @brief Assembles the stiffness matrix and system for a given potential field.
     * @param K Parallel bilinear form for the stiffness matrix.
     * @param gfc Coefficient for domain integration.
     * @param boundary Array of boundary markers.
     * @param potential Potential field grid function.
     * @param plf_B Linear form for the right-hand side.
     * @param matrix Resulting Hypre matrix.
     * @param hpv_X Solution vector.
     * @param hpv_B Right-hand-side vector.
     */
    void KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B);
    
    /**
     * @brief Configures and attaches a preconditioner to the conjugate gradient solver.
     * @param smoother Hypre preconditioner (e.g., Jacobi).
     * @param cg Conjugate gradient solver instance.
     * @param KMatrix Stiffness matrix to solve the system.
     */
    void PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix);
    
    /**
     * @brief Applies Dirichlet boundary conditions to the potential field.
     * @param dbc_Coef Coefficient for the boundary value.
     * @param Bv Boundary value.
     * @param phx Potential grid function to project the boundary conditions.
     * @param dbc_bdr Array of boundary markers.
     */
    void ImplementBoundaryConditions(mfem::ConstantCoefficient &dbc_Coef, double Bv, mfem::ParGridFunction &phx, mfem::Array<int> dbc_bdr);

    /**
     * @brief Creates a reaction field by scaling an input field with a factor.
     * @param Rx1 Input reaction field.
     * @param Rx2 Output reaction field (scaled version).
     * @param value Scaling factor.
     */
    void CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);
    
    /**
     * @brief Constructs the force term for the right-hand side of the system.
     * @param Rx2 Reaction field.
     * @param Fxx Linear form representing the force term.
     */
    void ForceTerm(mfem::ParGridFunction &Rx2, mfem::ParLinearForm &Fxx);
    
    /**
     * @brief Sets up the linear system considering boundary conditions.
     * @param K Stiffness matrix bilinear form.
     * @param boundary Array of boundary markers.
     * @param psx Input potential field.
     * @param plf_B Linear form for the right-hand side.
     * @param matrix Resulting Hypre matrix.
     * @param hpv_X Solution vector.
     * @param hpv_B Right-hand-side vector.
     * @param Coef Coefficient for boundary condition projection.
     * @param bdr Array of boundary attributes.
     */
    void ForceVector(mfem::ParBilinearForm &K, mfem::Array<int> boundary, mfem::ParGridFunction &psx, mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B, mfem::ConstantCoefficient &Coef, mfem::Array<int> &bdr);
    
    /**
     * @brief Calculates the error in the potential solution.
     * @param phx Potential field grid function.
     * @param cg_solver Conjugate gradient solver instance.
     * @param fterm Force term vector.
     * @param psx Auxiliary potential field.
     * @param error_X Local error accumulator.
     * @param gerror Global error (output).
     * @param gtPsx Total potential sum for normalization.
     */
    void ErrorCalculation(mfem::ParGridFunction &phx, mfem::CGSolver &cg_solver, mfem::HypreParVector &fterm, mfem::ParGridFunction &psx, double error_X, double &gerror, double gtPsx);

    double Vcell; ///< Cell voltage.

    int nE; ///< Number of elements in the mesh.
    int nC; ///< Number of nodes per element (corners).
    int nV; ///< Number of vertices in the mesh.

protected:
    
    mfem::ParMesh *pmesh; ///< Pointer to the parallel mesh.
    mfem::ParFiniteElementSpace *fespace; ///< Finite element space for the potentials.


private:

    mfem::ParGridFunction *Rxx; ///< Intermediate grid function for reactions.
    mfem::GridFunctionCoefficient *cXx; ///< Coefficient derived from grid functions.

    mfem::ParGridFunction *px0; ///< Grid function for potential values before iteration.
    mfem::HypreParVector X0; ///< Hypre vector for solving systems.

    mfem::ParGridFunction TmpF; ///< Temporary field for error calculations.

    const mfem::Vector& EVol; ///< Reference to element volumes from the MeshHandler.

};

#endif // POTENTIALS_HPP
