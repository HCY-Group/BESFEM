#ifndef CNP_HPP
#define CNP_HPP

#include "Concentrations_Base.hpp"

/**
 * @class CnP
 * @brief Derived class implementing particle concentration models for battery simulations.
 *
 * This class provides methods for initializing particle concentrations, 
 * performing time-stepping operations, and calculating reaction-related properties.
 */
class CnP : public Concentrations {

public:

    /**
     * @brief Constructor for the CnP class.
     * @param pmesh Pointer to the parallel mesh.
     * @param fespace Pointer to the finite element space.
     * @param mh Reference to the mesh handler.
     */
    CnP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);

    /**
     * @brief Initializes particle concentration values and solver components.
     * @param Cn Particle concentration grid function.
     * @param initial_value Initial concentration value.
     * @param psx Potential field grid function.
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Performs a single time step for particle concentration updates.
     * @param Rx Reaction grid function.
     * @param Cn Particle concentration grid function.
     * @param psx Potential field grid function.
     */
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    mfem::ParGridFunction *RxP; ///< Pointer to a grid function storing reaction values.



private:

    mfem::HypreParVector PsVc; ///< Vector for storing true degrees of freedom in the solid region.


    std::shared_ptr<mfem::CGSolver> Mp_solver; ///< Solver for the mass matrix.
    std::shared_ptr<mfem::HypreParMatrix> Mmatp; ///< Mass matrix for particle concentrations.
    mfem::HypreSmoother Mp_prec; ///< Preconditioner for the mass matrix solver.

    mfem::ParLinearForm ftPC; ///< Linear form for the force term related to particle concentration.

    std::shared_ptr<mfem::HypreParMatrix> Kmatp; ///< Stiffness matrix for diffusion calculations. 

    mfem::HypreParVector Fcb; ///< Right-hand-side vector for the system of equations.
    mfem::Array<int> boundary_dofs; ///< Array to store boundary degrees of freedom.
    mfem::HypreParVector X1v; ///< Temporary vector used during assembly.

    mfem::HypreParVector *CpV0; ///< Initial particle concentration values.
    mfem::HypreParVector *CpVn; ///< Particle concentration values at the next time step.
    mfem::HypreParVector *RHCp; ///< Right-hand-side vector at the current time step.
    mfem::HypreParMatrix *Tmatp; ///< System matrix for time-stepping.

};

#endif // CNP_HPP
