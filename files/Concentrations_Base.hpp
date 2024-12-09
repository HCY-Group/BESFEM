#ifndef CONCENTRATIONS_HPP
#define CONCENTRATIONS_HPP

// Public: Members can be accessed from anywhere. This is the default access modifier. 
// Protected: Members can be accessed within the class and by classes that inherit from that class. 
// Private: Members can only be accessed within the class that defines them.

#include "mfem.hpp"
#include "Mesh_Handler.hpp"

#include <memory>

/**
 * @class Concentrations
 * @brief Base class for managing concentration-related calculations in a finite element simulation.
 * 
 * This class provides a foundation for handling particle and electrolyte concentration
 * computations in a parallelized finite element framework. It includes functionality for:
 * - Initializing concentration values.
 * - Applying boundary conditions.
 * - Computing reactions, diffusivity, and related properties.
 */
class Concentrations {
public:

    /**
     * @brief Constructs a `Concentrations` object.
     * 
     * @param pmesh Pointer to the parallel mesh.
     * @param fespace Pointer to the parallel finite element space.
     * @param mh Reference to the mesh handler.
     */
    Concentrations(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    
    /**
     * @brief Virtual destructor for the `Concentrations` class.
     */
    virtual ~Concentrations() = default;

    /**
     * @brief Sets initial concentration values and prepares related parameters.
     * 
     * @param Cn Grid function for concentration to initialize.
     * @param initial_value Initial value for the concentration.
     * @param psx Grid function for psi values affecting the concentration.
     * @param perform_lithiation Flag to indicate whether to perform lithiation calculations.
     */
    void SetInitialValues(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);

    MeshHandler &mesh_handler; ///< Reference to the mesh handler object.
    mfem::ParGridFunction psi; ///< Grid function for psi values.                     
    mfem::ParGridFunction pse; ///< Grid function for pse values.

    mfem::ParGridFunction *CeT; ///< Pointer to the grid function for total concentration.

    // mfem::ParGridFunction *CnP;
    // mfem::ParGridFunction *CnE;

protected:
    mfem::ParMesh *pmesh; ///< Pointer to the parallel mesh.
    mfem::ParFiniteElementSpace *fespace; ///< Pointer to the finite element space.
    std::shared_ptr<mfem::HypreParMatrix> Mmat; ///< Pointer to the mass matrix.
    std::shared_ptr<mfem::CGSolver> solver; ///< Pointer to the conjugate gradient solver.
    mfem::HypreSmoother smoother; ///< Preconditioner for the solver.

    double infx; ///< Reaction current per unit width.
    double CeC = 0.0; ///< Local concentration sum.
    double gCeC = 0.0; ///< Global concentration sum.
    double CeAvg = 0.0; ///< Average concentration throughout the domain.
    double Ce0 = 0.001; ///< Initial reference concentration.

    /**
     * @brief Sets initial concentration values.
     * 
     * @param Cn Grid function for concentration.
     * @param initial_value Initial value to set.
     */
    void SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value);
    
    /**
     * @brief Sets up the solver for concentration calculations.
     * 
     * @param psx Grid function for psi values.
     * @param Mmat Matrix to store the mass matrix.
     * @param m_solver Conjugate gradient solver.
     * @param smoother Preconditioner for the solver.
     */
    void SetUpSolver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother);
    
    /**
     * @brief Applies Neumann boundary conditions.
     * 
     * @param psx Input grid function for psi values.
     * @param PGF Output grid function for boundary conditions.
     */
    void ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF);
    
    /**
     * @brief Creates a reaction term by scaling a grid function.
     * 
     * @param Rx1 Input reaction rate grid function.
     * @param Rx2 Output reaction rate grid function.
     * @param value Scaling factor for the reaction.
     */
    void CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);
    
    /**
     * @brief Computes a force term for the concentration system.
     * 
     * @param gfc Grid function for concentration.
     * @param Fxx Output linear form for force terms.
     * @param boundary Array specifying boundary conditions.
     * @param m Product coefficient for the force term.
     * @param apply_boundary_conditions Flag to apply boundary conditions.
     */
    void ForceTerm(mfem::ParGridFunction &gfc, mfem::ParLinearForm &Fxx, mfem::Array<int> boundary, mfem::ProductCoefficient m, bool apply_boundary_conditions);
    
    /**
     * @brief Computes the total reaction current.
     * 
     * @param Rx Grid function for reaction rate.
     * @param xCrnt Reference to the computed reaction current.
     */
    void TotalReaction(mfem::ParGridFunction &Rx, double xCrnt);
    
    /**
     * @brief Computes diffusivity values based on input parameters.
     * 
     * @param psx Grid function for psi values.
     * @param Cn Grid function for concentration values.
     * @param particle_electrolyte Flag to determine the computation equation.
     * @return Shared pointer to the computed grid function coefficient for diffusivity.
     */
    std::shared_ptr<mfem::GridFunctionCoefficient> Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte );
    
    /**
     * @brief Constructs the stiffness matrix for the concentration system.
     * 
     * @param boundary Array specifying boundary conditions.
     * @param Cn Grid function for concentration values.
     * @param Fxx Linear form for the force term.
     * @param Kmatx Shared pointer to the resulting stiffness matrix.
     * @param X1v Input vector for the linear system.
     * @param Fxb Right-hand side vector for the system.
     * @param cDx Coefficient representing diffusivity.
     */
    void KMatrix(mfem::Array<int> boundary, mfem::ParGridFunction &Cn, mfem::ParLinearForm &Fxx, std::shared_ptr<mfem::HypreParMatrix> &Kmatx, mfem::HypreParVector &X1v, mfem::HypreParVector &Fxb, mfem::GridFunctionCoefficient *cDx);
    
    /**
     * @brief Performs salt conservation calculations.
     * 
     * @param Cn Grid function for current concentration values.
     * @param psx Grid function for psi values.
     */
    void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
//     void Save(mfem::ParGridFunction &gf, const std::string &base_name);
    
    /**
     * @brief Performs lithiation calculations based on psi and concentration.
     * 
     * @param Cn Grid function for concentration values.
     * @param psx Grid function for psi values.
     */
    void LithiationCalculation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

private:
    // MeshHandler &mesh_handler;

    mfem::ParBilinearForm *M = nullptr; ///< Pointer to the bilinear form for mass matrix assembly.
    mfem::ParGridFunction *Ps_gf; ///< Pointer to the grid function for psi values.
    mfem::GridFunctionCoefficient *cP; ///< Pointer to the coefficient for psi grid function.

    mfem::Array<int> boundary_dofs; ///< Array of boundary degrees of freedom.

    mfem::ParGridFunction TmpF; ///< Temporary grid function for intermediate calculations.

    int nE; ///< Number of elements in the mesh.                                        
    int nC; ///< Number of corners per element.                                        
    int nV; ///< Number of vertices in the mesh.                                       

    const mfem::Vector& EVol; ///< Element volumes from MeshHandler.
    double gtPsi; ///< Total Psi from MeshHandler.
    double gtPse; ///< Total Pse from MeshHandler.

    mfem::ParGridFunction *Rxx; ///< Pointer to the grid function for reaction terms.
    mfem::GridFunctionCoefficient *cXx; ///< Pointer to the coefficient for reaction terms.

    // double infx;
    // double CeC = 0.0;
    // double gCeC = 0.0;
    // double CeAvg = 0.0;
    // double Ce0 = 0.001;

    double geCrnt;





};

#endif // CONCENTRATIONS_HPP
