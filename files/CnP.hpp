#ifndef CNP_HPP
#define CNP_HPP

#include "Concentrations_Base.hpp"

/**
 * @class CnP
 * @brief Manages the particle concentration calculations.
 * 
 * This class extends `Concentrations` to include specific functionality for:
 * - Initializing electrolyte concentrations.
 * - Time-stepping using Backward Euler.
 */
class CnP : public Concentrations {
public:
    /**
     * @brief Constructs a CnP object.
     * 
     * @param pm Pointer to the parallel mesh.
     * @param fe Pointer to the parallel finite element space.
     * @param mh Reference to the mesh handler.
     */
    CnP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    
    /**
     * @brief Initializes the particle concentration and solver setup.
     * 
     * @param Cn Grid function representing the initial concentration.
     * @param initial_value Initial value for the concentration.
     * @param psx Grid function for psi affecting concentration.
     * @param perform_lithiation Flag indicating whether to perform lithiation calculations.
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);

	 /**
     * @brief Performs a single time step for particle concentration calculations.
     * 
     * This method updates the electrolyte concentration based on reaction and diffusion
     * processes.
     * 
     * @param Rx Reaction rate grid function.
     * @param Cn Concentration grid function.
     * @param psx psi grid function.
     */
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    mfem::ParGridFunction *RxP; ///< Grid function for reaction rate.

private:

    mfem::HypreParVector PsVc;

    std::shared_ptr<mfem::CGSolver> Mp_solver; ///< Conjugate gradient solver for particle concentration.
    std::shared_ptr<mfem::HypreParMatrix> Mmatp; ///< Mass matrix for the particle.
    mfem::HypreSmoother Mp_prec; ///< Preconditioner for the CG solver.

    mfem::ParLinearForm ftPC; // force term particle concentration

    std::shared_ptr<mfem::HypreParMatrix> Kmatp; ///< Stiffness matrix for the particle.

    mfem::HypreParVector Fcb;///< Right-hand side vector for force terms.
    mfem::Array<int> boundary_dofs; ///< Array of boundary degrees of freedom.
    mfem::HypreParVector X1v;

    mfem::HypreParVector *CpV0; ///< Initial particle concentration vector.
    mfem::HypreParVector *CpVn; ///< Updated particle concentration vector.
    mfem::HypreParVector *RHCp; ///< Right-hand side vector for Backward Euler updates.
    mfem::HypreParMatrix *Tmatp; ///< Matrix for Backward Euler.

};

#endif // CNP_HPP
