#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"
#include "Constants.hpp"
#include "SimTypes.hpp"

/**
 * @file PotE.hpp
 * @brief Electrolyte potential solver for BESFEM.
 *
 * Implements initialization, conductivity/diffusivity assembly, and 
 * time-stepping for the electrolyte potential field φ_E. This class 
 * extends the generic PotentialBase interface and incorporates 
 * electrolyte-specific boundary conditions, phase-field masking, and 
 * reaction coupling.
 */

/**
 * @brief Dirichlet boundary value for the electrolyte potential φ_E.
 */
extern double BvE;

/**
 * @class PotE
 * @brief Electrolyte-phase potential solver.
 *
 * Responsibilities:
 * - Initialize φ_E using user-specified initial conditions
 * - Assemble conductivity and diffusivity operators (κ_e, D_e)
 * - Assemble reaction/forcing terms
 * - Solve the linear system each timestep (Hypre CG + AMG)
 * - Provide optional coupling with concentration true-DoF vectors
 * - Report global error metrics for convergence monitoring
 */
class PotE : public PotentialBase {

public:

    /**
     * @brief Construct the electrolyte potential solver.
     *
     * Stores references to geometry, domain parameters, and boundary 
     * conditions and initializes FEM helper utilities.
     *
     * @param geo   Geometry/mesh handler.
     * @param para  Domain parameter container.
     * @param bc    Boundary condition handler.
     * @param mode  Cell mode (HALF or FULL).
     */
    PotE(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::CellMode mode);

    Initialize_Geometry  &geometry;            ///< Geometry and mesh infrastructure.
    Domain_Parameters    &domain_parameters;   ///< Material and phase-field parameters.
    BoundaryConditions   &boundary_conditions; ///< Electrode/BC configuration.
    sim::CellMode         mode_;               ///< Half/full-cell selection mode.

    FEMOperators fem; ///< FEM operator assembly utilities.
    Utils utils;      ///< Miscellaneous helpers (clamping, reductions, etc.).

    double BvE;       ///< Dirichlet boundary value for φ_E (if applicable).

    /**
     * @brief Initialize electrolyte-phase potential φ_E.
     *
     * Tasks performed:
     * - Seed initial potential field.
     * - Build conductivity/diffusivity operators.
     * - Setup solver and Hypre AMG preconditioner.
     *
     * @param Cn            Electrolyte concentration field.
     * @param initial_value Initial scalar value for φ_E.
     * @param psx           Phase-field mask ψ_E.
     */
    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Assemble the linear system for φ_E.
     *
     * This overload includes a true-DoF concentration vector CeVn
     * which is used in high-fidelity coupling models (e.g., electrolyte transport).
     *
     * @param Cn        Electrolyte concentration field.
     * @param psx       Phase-field mask ψ_E.
     * @param potential Potential field φ_E (in/out).
     * @param CeVn      True-DoF concentration vector for coupling.
     */
    void AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                        mfem::ParGridFunction &potential, mfem::HypreParVector  &CeVn);

    /**
     * @brief Assemble the linear system for φ_E.
     *
     * @param Cn        Electrolyte concentration field.
     * @param psx       Phase-field mask ψ_E.
     * @param potential Potential field φ_E (in/out).
     */
    void AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                        mfem::ParGridFunction &potential);

    /**
     * @brief Advance φ_E by one timestep.
     *
     * Recomputes forcing terms, solves for the updated potential, and computes 
     * a global MPI-reduced error metric.
     *
     * @param Rx     Reaction source field.
     * @param phx    Potential field φ_E (in/out).
     * @param psx    Phase mask ψ_E.
     * @param gerror Output: global RMS/L2 error value.
     */
    void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx,
                         mfem::ParGridFunction &psx, double &gerror);

    /**
     * @brief Advance φ_E using two reaction fields (multi-electrode/full-cell cases).
     *
     * @param Rx1    Reaction source term #1.
     * @param Rx2    Reaction source term #2.
     * @param phx    Potential field φ_E (in/out).
     * @param psx    Phase mask ψ_E.
     * @param gerror Output: global RMS/L2 error value.
     */
    void UpdatePotential(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2,
                         mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);

private:

    // -------------------------------------------------------------------------
    // Spaces / BC handling
    // -------------------------------------------------------------------------
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space for φ_E.

    mfem::Array<int> dbc_bdr;           ///< Dirichlet boundary attribute markers.
    mfem::Array<int> ess_tdof_list_potE;///< Essential true DOFs for φ_E BCs.

    // -------------------------------------------------------------------------
    // Physical constants (transport coefficients)
    // -------------------------------------------------------------------------
    double tc1 = (2 * Constants::t_minus - 1.0) / (2 * Constants::t_minus * (1.0 - Constants::t_minus));
    double tc2 = Constants::Cst1 / (2 * Constants::t_minus * (1.0 - Constants::t_minus));

    double dffe  = 0.0; ///< Temporary scratch value.
    double gtPse = 0.0; ///< Global total ψ_E.

    // -------------------------------------------------------------------------
    // Solver and preconditioner
    // -------------------------------------------------------------------------
    mfem::CGSolver                       cgPE_solver; ///< CG solver for φ_E.
    std::unique_ptr<mfem::HypreBoomerAMG> Mpe;         ///< AMG preconditioner.

    // -------------------------------------------------------------------------
    // Operators (bilinear/linear forms and matrices)
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParBilinearForm> Kl1; ///< Conductivity form.
    std::unique_ptr<mfem::ParBilinearForm> Kl2; ///< Diffusivity / extra form.
    std::unique_ptr<mfem::ParLinearForm>    Bl2; ///< Reaction/forcing form.

    mfem::HypreParMatrix Kml; ///< Conductivity matrix.
    mfem::HypreParMatrix Kdm; ///< Diffusivity matrix.

    // -------------------------------------------------------------------------
    // Coefficients and material fields
    // -------------------------------------------------------------------------
    mfem::ParGridFunction kpl; ///< Conductivity κ_e(x).
    mfem::ParGridFunction Dmp; ///< Diffusivity D_e(x).

    mfem::GridFunctionCoefficient cKe; ///< κ_e wrapper.
    mfem::GridFunctionCoefficient cDm; ///< D_e wrapper.
    mfem::GridFunctionCoefficient cRe; ///< Reaction coefficient wrapper.

    mfem::ParGridFunction RpE; ///< Reaction field (grid function).
    mfem::ParGridFunction pE0; ///< Previous φ_E (grid scratch).

    // -------------------------------------------------------------------------
    // Vectors and RHS data
    // -------------------------------------------------------------------------
    mfem::ParLinearForm  B1t; ///< RHS (domain + boundary forcing).
    mfem::ParLinearForm Flt;   ///< Force term linear form.
    mfem::HypreParVector X1v; ///< Solution vector.
    mfem::HypreParVector B1v; ///< RHS vector.
    mfem::HypreParVector Flb; ///< Additional RHS buffer.
    mfem::HypreParVector Xe0; ///< Scratch solution vector.
    mfem::HypreParVector RHSl;///< RHS vector (assembled).
    mfem::HypreParVector LpCe;///< Concentration true-DoF vector (coupling).
    mfem::HypreParVector CeVn; ///< Concentration at next time step



    mfem::Array<int> boundary_dofs; ///< Boundary DOF list.

    // -------------------------------------------------------------------------
    // Anchor/pinning (to fix gauge freedom)
    // -------------------------------------------------------------------------
    mfem::ParGridFunction phE_bc; ///< BC/pinning helper field.
    bool anchor_set = false;       ///< Whether a gauge anchor has been chosen.
    int  myid       = 0;           ///< MPI rank.
    mfem::Array<int> ess_tdof_potE;

    bool pin  = false; ///< Whether a pin exists.
    int  rkpp = -1;    ///< Rank owning the pinned DOF.

    mfem::Array<int> dbc_e_bdr; ///< Additional Dirichlet markers (east).

    /**
     * @brief Recompute electrolyte conductivity κ_e and diffusivity D_e.
     *
     * Applies ψ_E masking and concentration dependence.
     *
     * @param Cn  Electrolyte concentration field.
     * @param psx Phase-field ψ_E.
     */
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
};

#endif // POTE_HPP
