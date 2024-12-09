#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"


/**
 * @class CnE
 * @brief Manages the electrolyte concentration calculations.
 * 
 * This class extends `Concentrations` to include specific functionality for:
 * - Initializing electrolyte concentrations.
 * - Time-stepping using Crank-Nicolson.
 * - Applying boundary conditions.
 */
class CnE : public Concentrations {
public:
    /**
     * @brief Constructs a CnE object.
     * 
     * @param pm Pointer to the parallel mesh.
     * @param fe Pointer to the parallel finite element space.
     * @param mh Reference to the mesh handler.
     */
    CnE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    
    /**
     * @brief Initializes the electrolyte concentration and solver setup.
     * 
     * @param Cn Grid function representing the initial concentration.
     * @param initial_value Initial value for the concentration.
     * @param psx Grid function for psi affecting concentration.
     * @param perform_lithiation Flag indicating whether to perform lithiation calculations.
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);


	 /**
     * @brief Performs a single time step for electrolyte concentration calculations.
     * 
     * This method updates the electrolyte concentration based on reaction and diffusion
     * processes, including boundary condition application.
     * 
     * @param Rx Reaction rate grid function.
     * @param Cn Concentration grid function.
     * @param psx psi grid function.
     */
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    
    mfem::ParGridFunction *RxE; ///< Grid function for reaction rate.
    mfem::ParLinearForm ftE; ///< Linear form for reaction force terms.
    mfem::Array<int> nbc_w_bdr; ///< West Side Neumann Boundary Conditions.
    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; 
    mfem::Array<int> boundary_dofs; ///< Array of boundary degrees of freedom.

    static mfem::HypreParVector *CeVn; ///< Updated electrolyte concentration vector.
    static mfem::HypreParVector* GetCeVn() { return CeVn; } // static variable to be used in reaction



private:


    mfem::ParGridFunction *PeR; ///< Grid function for reaction potential.

    std::shared_ptr<mfem::CGSolver> Me_solver; ///< Conjugate gradient solver for electrolyte concentration.
    std::shared_ptr<mfem::HypreParMatrix> Mmate; ///< Mass matrix for the electrolyte.
    mfem::HypreSmoother Me_prec; ///< Preconditioner for the CG solver.

    double eCrnt; ///< Current reaction rate.

    std::shared_ptr<mfem::HypreParMatrix> Kmate; ///< Stiffness matrix for the electrolyte.

    mfem::HypreParVector Feb; ///< Right-hand side vector for force terms.
    mfem::HypreParVector X1v;

    mfem::HypreParVector *CeV0; ///< Initial electrolyte concentration vector.
    mfem::HypreParVector *RHCe; ///< Right-hand side vector for Crank-Nicolson updates.
    mfem::HypreParMatrix *TmatR; ///< Right-hand side Crank-Nicolson matrix.
    mfem::HypreParMatrix *TmatL; ///< Left-hand side Crank-Nicolson matrix.



};

#endif // CNE_HPP
