#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"
#include "../inputs/Constants.hpp"
#include "SimTypes.hpp"

/**
 * @file PotE.hpp
 * @brief Electrolyte potential (PotE) solver for battery simulations.
 *
 * Extends the Potentials base class to assemble/solve the electrolyte
 * potential system, including conductivity/diffusivity coupling,
 * time-stepping helpers, and error reporting.
 */

/**
 * @brief Boundary value for the electrolyte potential.
 */
extern double BvE;

/**
 * @defgroup potentials Potential Modules
 * @brief Classes that assemble and advance potential fields.
 * @{
 */

/**
 * @class PotE
 * @brief Electrolyte potential solver (electrolyte phase).
 * @ingroup potentials
 *
 * @details
 * - Builds conductivity/diffusivity operators in MFEM/Hypre.
 * - Handles Dirichlet boundary conditions on marked boundaries.
 * - Provides time-stepping helpers and residual/error evaluation.
 */
class PotE : public Potentials {

public:

    /**
     * @brief Construct the electrolyte potential solver.
     * @param geo  Geometry/space container (mesh, FESpaces, BC markers).
     * @param para Domain/physics parameters (constants, totals, dt).
     */
    PotE(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc);

    Initialize_Geometry &geometry; ///< Reference to geometry initialization
    Domain_Parameters &domain_parameters; ///< Reference to domain parameters
    BoundaryConditions &boundary_conditions;


    double BvE; ///< Boundary value for electrolyte potential

    /**
     * @brief Initialize electrolyte potential operators and state.
     *
     * Assembles conductivity/diffusivity forms, configures solver/preconditioner,
     * and seeds the potential field.
     *
     * @param Cn            Electrolyte concentration field (input).
     * @param initial_value Initial scalar value for the potential field.
     * @param psx           Phase mask ψ_E used for weighting/BCs.
     */
    void Initialize(sim::CellMode mode, mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Advance electrolyte potential for one step (assemble/solve).
     *
     * Uses current concentration and phase mask to update operators and solve
     * for the electrolyte potential.
     *
     * @param Cn        Electrolyte concentration field (input).
     * @param psx       Phase mask ψ_E (input).
     * @param potential Electrolyte potential grid function (in/out).
     * @param CeVn      True-DoF concentration vector (input), for coupling.
     */
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn);

    /**
     * @brief Advance electrolyte potential for one step (assemble/solve).
     *
     * Uses current concentration and phase mask to update operators and solve
     * for the electrolyte potential.
     *
     * @param Cn        Electrolyte concentration field (input).
     * @param psx       Phase mask ψ_E (input).
     * @param potential Electrolyte potential grid function (in/out).
     */
    void TimeStep(sim::CellMode mode, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential);
    
    
    
    /**
     * @brief One higher-level advance (e.g., residual/error evaluation).
     *
     * @param Rx     Reaction source field (input).
     * @param phx    Electrolyte potential field (in/out).
     * @param psx    Phase mask ψ_E (input).
     * @param gerror Global error/residual (output).
     */
    void Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);


    /**
     * @brief One higher-level advance with reaction field update.
     * @param Rx1    Input reaction source field (input).
     * @param Rx2    Input reaction source field (input).
     * @param phx    Electrolyte potential field (in/out).
     * @param psx    Phase mask ψ_E (input).
     * @param gerror Global error/residual (output).
     */
    void Advance(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);
   

private:

    // ---------- Spaces / BC bookkeeping ----------
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.
    mfem::Array<int> dbc_bdr;       ///< Dirichlet boundary markers (west).
    mfem::Array<int> ess_tdof_list_potE; ///< Essential true DOFs on west boundary.

    // ---------- Physical constants (derived combos) ----------
    double tc1 = (2*Constants::t_minus - 1.0) / (2*Constants::t_minus*(1.0 - Constants::t_minus)); ///< Transport factor 1.
    double tc2 = Constants::Cst1 / (2*Constants::t_minus*(1.0 - Constants::t_minus));               ///< Transport factor 2.

    double dffe = 0.0; ///< Scratch for conductivity/reaction calculations.
    double gtPse = 0.0; ///< Global total for Pse (electrolyte phase).

    // ---------- Solver / preconditioner ----------
    mfem::CGSolver       cgPE_solver; ///< Conjugate gradient solver.
    std::unique_ptr<mfem::HypreBoomerAMG> Mpe; ///< BoomerAMG preconditioner (built once with KmP)


    // ---------- Forms / matrices ----------
    std::unique_ptr<mfem::ParBilinearForm> Kl1; ///< Conductivity bilinear form.
    std::unique_ptr<mfem::ParBilinearForm> Kl2; ///< Additional conductivity/diffusivity form.
    std::unique_ptr<mfem::ParLinearForm>    Bl2; ///< Reaction/load linear form.

    mfem::HypreParMatrix Kml; ///< Stiffness matrix (conductivity).
    mfem::HypreParMatrix Kdm; ///< Stiffness matrix (diffusivity).
    mfem::HypreParVector CeVn; ///< Concentration at next time step


    // ---------- Coefficients / fields ----------
    mfem::ParGridFunction kpl; ///< Conductivity field κ_e(x).
    mfem::ParGridFunction Dmp; ///< Diffusivity field D_e(x).
    mfem::GridFunctionCoefficient cKe; ///< Coefficient wrapping κ_e.
    mfem::GridFunctionCoefficient cDm; ///< Coefficient wrapping D_e.
    mfem::GridFunctionCoefficient cRe; ///< Coefficient wrapping reaction term.
    mfem::ParGridFunction RpE; ///< Reaction field for electrolyte potential.

    // ---------- Vectors / RHS ----------
    mfem::ParLinearForm  B1t; ///< RHS linear form (domain/boundary forcing).
    mfem::HypreParVector X1v; ///< Solution vector (true DOFs).
    mfem::HypreParVector B1v; ///< RHS vector.
    mfem::HypreParVector Flb; ///< Auxiliary RHS vector.
    mfem::HypreParVector Xe0; ///< Scratch solution vector.
    mfem::ParGridFunction pE0; ///< Electrolyte potential (grid) scratch.
    mfem::HypreParVector RHSl; ///< RHS vector (assembled).

    mfem::ParLinearForm Flt;   ///< Force term linear form.
    mfem::Array<int>     boundary_dofs; ///< Boundary DOFs list.

    mfem::HypreParVector LpCe; ///< Concentration true-DoF vector (coupling).

    mfem::ParGridFunction phE_bc;         // reused BC vector for gauge pin
    bool anchor_set = false;               // track if anchor has been set
    int myid;
    mfem::Array<int> ess_tdof_potE;

    bool pin = false;
    int rkpp = -1;
    mfem::Array<int> dbc_e_bdr; ///< Array marking Dirichlet boundary attributes


    

    /**
     * @brief Recompute electrolyte conductivity/diffusivity fields and coefficients.
     * @param Cn  Electrolyte concentration field.
     * @param psx Phase mask ψ_E.
     */
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

};

/** @} */ // end group potentials


#endif // POTE_HPP