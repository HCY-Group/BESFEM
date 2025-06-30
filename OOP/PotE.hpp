#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"
#include "../code/Constants.hpp"
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

public:

    PotE(Initialize_Geometry &geo, Domain_Parameters &para);

    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;

    double BvE; ///< Boundary value for electrolyte potential


    // /**
    //  * @brief Initializes the electrolyte potential field and solver
    //  * @param ph Electrolyte potential field
    //  * @param initial_value Initial potential value
    //  */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    // /**
    //  * @brief Performs time-stepping for the electrolyte potential
    //  * @param Cn Electrolyte concentration field
    //  * @param psx Psi potential field
    //  * @param phx Electrolyte potential field
    //  */
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn);
    
    // /**
    //  * @brief Calculates the global error in the electrolyte potential solution
    //  * @param Rx Reaction field
    //  * @param phx Electrolyte potential field
    //  * @param psx Psi potential field
    //  * @param gerror Global error (output)
    //  */
    // void CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);

    // std::unique_ptr<mfem::CGSolver> cgPE_solver; ///< Variable for the conjugate gradient solver

    // mfem::Array<int> ess_tdof_list_w; ///< List of essential true degrees of freedom for Dirichlet boundary conditions
    // mfem::ConstantCoefficient dbc_w_Coef; ///< Coefficient for Dirichlet boundary conditions
    // mfem::Array<int> dbc_w_bdr; ///< Array marking Dirichlet boundary attributes

    // double error_E = 1.0; ///< Local error in the electrolyte potential solution


private:

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space

    mfem::Array<int> dbc_w_bdr; ///< Array marking Dirichlet boundary attributes
    mfem::Array<int> ess_tdof_list_w; ///< List of essential true degrees of freedom for Dirichlet boundary conditions

    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus)); ///< Constant for conductivity and diffusivity calculations
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1; ///< Constant for conductivity and diffusivity calculations

    double dffe; ///< Temporary variable for conductivity and reaction calculations
    
    std::shared_ptr<mfem::CGSolver> cgPE_solver;

    mfem::ParLinearForm B1t; ///< Linear form for the right-hand side
    mfem::HypreParVector X1v; ///< Solution vector
    mfem::HypreParVector B1v; ///< Right-hand-side vector
    mfem::HypreParVector Flb; ///< Right-hand-side vector
    mfem::HypreSmoother Mpe; ///< Preconditioner for the solver

    double gtPse; ///< Total Psi from MeshHandler

    mfem::ParGridFunction kpl; ///< Conductivity field for electrolyte potential
    mfem::ParGridFunction Dmp; ///< Diffusivity field for electrolyte potential
    mfem::GridFunctionCoefficient cKe; ///< Coefficient for the conductivity field
    mfem::GridFunctionCoefficient cDm; ///< Coefficient for the diffusivity field
    mfem::GridFunctionCoefficient cRe; ///< Coefficient for the reaction field
    mfem::ParGridFunction RpE; ///< Reaction field for the electrolyte potential

    std::unique_ptr<mfem::ParBilinearForm> Kl1; ///< Bilinear form for electrolyte potential conductivity
    std::unique_ptr<mfem::ParBilinearForm> Kl2; ///< Bilinear form for electrolyte potential conductivity
    std::shared_ptr<mfem::HypreParMatrix> Kml; ///< Stiffness matrix for electrolyte potential conductivity
    std::shared_ptr<mfem::HypreParMatrix> Kdm; /// Stiffness matrix for diffusivity
    std::unique_ptr<mfem::ParLinearForm> Bl2; ///< Linear form for the reaction term

    mfem::ParLinearForm Flt; ///< Linear form for the force term in electrolyte potential calculations
    mfem::Array<int> boundary_dofs; ///< Array of boundary degrees of freedom   

    mfem::HypreParVector LpCe; ///< Pointer to the vector for concentration degrees of freedom

    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx); ///< Computes electrolyte conductivity and diffusivity

};

#endif // POTE_HPP