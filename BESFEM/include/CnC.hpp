#ifndef CNC_HPP
#define CNC_HPP

#include "Concentrations_Base.hpp"

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class CnC
 * @brief Derived class implementing particle concentration models for battery simulations.
 *
 * This class provides methods for initializing particle concentrations, 
 * performing time-stepping operations, and calculating reaction-related properties.
 */
class CnC : public Concentrations {

public:

    /**
     * @brief Constructor for the CnC class
     * @param pm Pointer to the parallel mesh
     * @param fe Pointer to the finite element space
     * @param mh Reference to the mesh handler
     */
    CnC(Initialize_Geometry &geo, Domain_Parameters &para);


    /**
     * @brief Initializes particle concentration values and solver components
     * @param Cn Particle concentration grid function
     * @param initial_value Initial concentration value
     * @param psx Potential field grid function
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Performs a single time step for particle concentration updates
     * @param Rx Reaction grid function
     * @param Cn Particle concentration grid function
     * @param psx Potential field grid function
     */
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);



private:

    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;

    mfem::HypreParVector PsVc; ///< Vector for storing true degrees of freedom in the solid region
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space

    mfem::ParGridFunction Dp; ///< Grid function for particle diffusivity
    mfem::ParGridFunction RxP; ///< Pointer to a grid function storing reaction values

    mfem::CGSolver Mp_solver; ///< Solver for the mass matrix
    std::shared_ptr<mfem::HypreParMatrix> Mmatp; ///< Mass matrix for particle concentrations
    mfem::HypreSmoother Mp_prec; ///< Preconditioner for the mass matrix solver

    mfem::ParLinearForm ftPC; ///< Linear form for the force term related to particle concentration

    std::shared_ptr<mfem::HypreParMatrix> Kmatp; ///< Stiffness matrix for diffusion calculations 

    mfem::HypreParVector Fcb; ///< Vector for storing the force term contributions
    mfem::ParLinearForm Fct; ///< Linear form for the force term related to particle concentrations

    mfem::Array<int> boundary_dofs; ///< Array to store boundary degrees of freedom
    mfem::HypreParVector X1v; ///< Temporary vector used during assembly

    std::shared_ptr<mfem::HypreParVector> CpV0; ///< Initial particle concentration values
    std::shared_ptr<mfem::HypreParVector> CpVn; ///< Particle concentration values at the next time step
    std::shared_ptr<mfem::HypreParVector> RHCp; ///< Right-hand-side vector at the current time step
    std::unique_ptr<mfem::HypreParMatrix> Tmatp; ///< System matrix for time-stepping

    std::unique_ptr<mfem::ParBilinearForm> Mt; ///< Mass matrix for particle concentrations

    mfem::GridFunctionCoefficient cAp; ///< Coefficient for the reaction term
    mfem::GridFunctionCoefficient cDp; ///< Coefficient for the diffusivity term
    std::unique_ptr<mfem::ParLinearForm> Bc2; ///< Linear form for the force term related to particle concentrations
    std::unique_ptr<mfem::ParBilinearForm> Kc2; ///< Stiffness form for particle potential
};

#endif // CNC_HPP
