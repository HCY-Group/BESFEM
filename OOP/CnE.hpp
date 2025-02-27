#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"

/**
 * @class CnE
 * @brief Derived class implementing electrolyte concentration models for battery simulations.
 *
 * This class provides methods for initializing electrolyte concentrations, 
 * performing time-stepping operations, and managing reaction and boundary conditions.
 */

class CnE : public Concentrations {
public:

    /**
     * @brief Constructor for the CnE class
     * @param pm Pointer to the parallel mesh
     * @param fe Pointer to the finite element space
     * @param mh Reference to the mesh handler
     */
    CnE(Initialize_Geometry &geo, Domain_Parameters &para);

    // virtual ~CnE();

    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;

    /**
     * @brief Initializes electrolyte concentration values and solver components
     * @param Cn Electrolyte concentration grid function
     * @param initial_value Initial concentration value
     * @param psx Potential field grid function
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Performs a single time step for electrolyte concentration updates
     * @param Rx Reaction grid function
     * @param Cn Electrolyte concentration grid function
     * @param psx Potential field grid function
     */
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    // mfem::ParGridFunction *RxE; ///< Pointer to a grid function storing reaction values for the electrolyte
    std::unique_ptr<mfem::ParGridFunction> RxE;

    mfem::ParLinearForm ftE; ///< Linear form for the force term related to electrolyte concentration
    mfem::Array<int> nbc_w_bdr; ///< Boundary markers for Neumann boundary conditions
    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; ///< Product coefficient for Neumann boundary conditions
    mfem::Array<int> boundary_dofs; ///< Array to store boundary degrees of freedom

    // static mfem::HypreParVector *CeVn; ///< Static variable to store the electrolyte concentration vector at the next time step
    std::shared_ptr<mfem::HypreParVector> CeVn;

    // /**
    //  * @brief Getter for the static variable CeVn
    //  * @return Pointer to the electrolyte concentration vector at the next time step
    //  */
    // static mfem::HypreParVector* GetCeVn() { return CeVn; } // static variable to be used in reaction

private:

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space

    // mfem::ParGridFunction *PeR; ///< Pointer to a grid function storing reaction potential values
    std::unique_ptr<mfem::ParGridFunction> PeR;



    std::shared_ptr<mfem::CGSolver> Me_solver; ///< Solver for the mass matrix
    std::shared_ptr<mfem::HypreParMatrix> Mmate; ///< Mass matrix for electrolyte concentrations
    mfem::HypreSmoother Me_prec; ///< Preconditioner for the mass matrix solver

    double eCrnt; ///< Current reaction value for the electrolyte

    std::shared_ptr<mfem::HypreParMatrix> Kmate; ///< Stiffness matrix for diffusion calculations

    mfem::HypreParVector Feb; ///< Right-hand-side vector for the system of equations
    mfem::HypreParVector X1v; ///< Temporary vector used during assembly

    // mfem::HypreParVector *CeV0; ///< Initial electrolyte concentration values
    // mfem::HypreParVector *RHCe; ///< Right-hand-side vector at the current time step
    // mfem::HypreParMatrix *TmatR; ///< System matrix for the right-hand-side calculation (Crank-Nicolson)
    // mfem::HypreParMatrix *TmatL; ///< System matrix for the left-hand-side calculation (Crank-Nicolson)

    std::unique_ptr<mfem::HypreParMatrix> TmatR;
    std::unique_ptr<mfem::HypreParMatrix> TmatL;

    std::shared_ptr<mfem::HypreParVector> CeV0;
    std::shared_ptr<mfem::HypreParVector> RHCe;


    std::shared_ptr<mfem::ParBilinearForm> eKx2;



};

#endif // CNE_HPP
