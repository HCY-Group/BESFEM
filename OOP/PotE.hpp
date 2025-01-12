#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"

/**
 * @brief Boundary value for the electrolyte potential.
 */
extern double BvE;

/**
 * @class PotE
 * @brief Derived class implementing the electrolyte potential model for battery simulations.
 *
 * This class extends the `Potentials` base class to handle the computation
 * of electrolyte potentials, including initialization, time-stepping, and error calculations.
 */
class PotE : public Potentials {

    std::unique_ptr<mfem::ParBilinearForm> Kl1; ///< Unique pointer for bilinear forms used in matrix assembly
    std::unique_ptr<mfem::ParBilinearForm> Kl2; ///< Unique pointer for bilinear forms used in matrix assembly



public:

    /**
     * @brief Constructor for the PotE class
     * @param pm Pointer to the parallel mesh
     * @param fe Pointer to the finite element space
     * @param mh Reference to the mesh handler
     */
    PotE(mfem::ParMesh* pm, mfem::ParFiniteElementSpace* fe, MeshHandler &mh);

    /**
     * @brief Initializes the electrolyte potential field and solver
     * @param ph Electrolyte potential field
     * @param initial_value Initial potential value
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);

    /**
     * @brief Performs time-stepping for the electrolyte potential
     * @param Cn Electrolyte concentration field
     * @param psx Psi potential field
     * @param phx Electrolyte potential field
     */
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx);
    
    /**
     * @brief Calculates the global error in the electrolyte potential solution
     * @param Rx Reaction field
     * @param phx Electrolyte potential field
     * @param psx Psi potential field
     * @param gerror Global error (output)
     */
    void CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);

    static mfem::CGSolver *cgPE_solver; ///< Static variable for the conjugate gradient solver
    
    /**
     * @brief Getter for the static conjugate gradient solver instance
     * @return Pointer to the conjugate gradient solver
     */
    static mfem::CGSolver *GetcgPEsolver() { return cgPE_solver; }

    mfem::Array<int> ess_tdof_list_w; ///< List of essential true degrees of freedom for Dirichlet boundary conditions
    mfem::ConstantCoefficient dbc_w_Coef; ///< Coefficient for Dirichlet boundary conditions
    mfem::Array<int> dbc_w_bdr; ///< Array marking Dirichlet boundary attributes

    double error_E = 1.0; ///< Local error in the electrolyte potential solution


private:

    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus)); ///< Constant for conductivity and diffusivity calculations
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1; ///< Constant for conductivity and diffusivity calculations

    double dffe; ///< Temporary variable for conductivity and reaction calculations
    
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx); ///< Computes electrolyte conductivity and diffusivity

    mfem::ParGridFunction *Dmp; ///< Diffusivity field - D minus plus
    mfem::ParGridFunction *kpl; ///< Conductivity field

    mfem::HypreParMatrix Kdm; ///< Stiffness matrix for diffusivity
    mfem::HypreParMatrix KmE; ///< Stiffness matrix for conductivity

    mfem::ParLinearForm B1t; ///< Linear form for the right-hand side
    mfem::HypreParVector X1v; ///< Solution vector
    mfem::HypreParVector B1v; ///< Right-hand-side vector
    mfem::HypreParVector RHSl; ///< Residual vector

    mfem::ParFiniteElementSpace *fespace; // Declare the member variable fe

    mfem::HypreParVector *CeVn; ///< Concentration at the current time step
    mfem::HypreParVector *LpCe; ///< Diffusion-related vector

    mfem::HypreSmoother Mpe; ///< Preconditioner for the solver

    Array<int> boundary_dofs; ///< Array of boundary degrees of freedom

    mfem::GridFunctionCoefficient cKp; ///< Coefficient for the stiffness matrix

    mfem::ParGridFunction *RpE; ///< Reaction field for the electrolyte
    mfem::ParLinearForm ftPotE; // force term particle electrolyte
    mfem::HypreParVector Flb; ///< Vector for assembled force term
    
    double gtPse; ///< Total potential sum from MeshHandler

};

#endif // POTE_HPP