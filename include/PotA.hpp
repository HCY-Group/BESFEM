#ifndef POTA_HPP
#define POTA_HPP

/**
 * @file PotA.hpp
 * @brief Anode potential solver for BESFEM.
 *
 * Implements initialization, conductivity assembly, and time-stepping for
 * the solid active-material (particle) potential field φ_A. This class
 * extends the generic PotentialBase interface to handle anode-specific
 * boundary conditions, conductivity fields, and reaction coupling.
 */

#include "Potentials_Base.hpp"

/**
 * @brief Boundary value for the anode (particle) potential.
 *
 */
extern double BvA;

/**
 * @class PotA
 * @brief Solid-phase anode potential solver.
 *
 * Responsibilities:
 * - Initialize φ_A and apply Dirichlet/Neumann boundary conditions
 * - Assemble conductivity-weighted stiffness matrices
 * - Assemble reaction/forcing terms
 * - Solve the linear system each timestep (Hypre CG + AMG)
 * - Compute and return global error metrics
 */
class PotA : public PotentialBase {

public:

    /**
     * @brief Construct the anode potential solver.
     *
     * Stores references to geometry, domain parameters, and boundary
     * conditions and creates helper classes for FEM assembly.
     *
     * @param geo  Geometry/mesh handler (serial + parallel).
     * @param para Domain parameter container.
     * @param bc   Boundary condition handler.
     */
    PotA(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc);

    Initialize_Geometry &geometry;          ///< Geometry and mesh infrastructure.
    Domain_Parameters   &domain_parameters; ///< Material and phase-field parameters.
    BoundaryConditions  &boundary_conditions; ///< Electrode/BC configuration.

    FEMOperators fem; ///< FEM operator assembly utilities.
    Utils utils;      ///< Utility helpers (clamping, reduction, etc.).

    double BvA; ///< Dirichlet boundary value for φ_A (if applicable).

    /**
     * @brief Initialize solid-phase potential φ_A.
     *
     * Performs:
     * - Initial assignment of φ_A
     * - Setup of essential boundary DOFs
     * - Assembly of conductivity operators
     * - Preconditioning setup
     *
     * @param ph            Potential field to initialize.
     * @param initial_value Scalar initial potential value.
     * @param psx           Phase-field mask ψ_A.
     */
    void SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx) override;

    /**
     * @brief Assemble linear system for the anode potential.
     *
     * Builds:
     * - conductivity-weighted stiffness matrix,
     * - reaction linear form,
     * - boundary condition modifications.
     *
     * @param Cn        Solid concentration field.
     * @param psx       Phase-field mask ψ_A.
     * @param potential φ_A grid function (in/out).
     */
    void AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential) override;

    /**
     * @brief Advance φ_A by one timestep.
     * 
     * Updates:
     * - reaction field
     * - forcing terms
     * - global error `gerror`
     *
     * @param Rx     Reaction source term.
     * @param phx    Potential field φ_A (in/out).
     * @param psx    Phase-field ψ_A.
     * @param gerror Global RMS/L2 error metric (MPI reduced).
     */
    void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx,
                         double &gerror) override;

private:

    /**
     * @brief Compute solid-phase conductivity κ_p(C, ψ_A).
     *
     * Uses concentration-dependent or region-dependent conductivity rules
     * to populate the κ_p grid function.
     *
     * @param Cn  Concentration field.
     * @param psx Phase-field ψ_A.
     */
    void ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space for φ_A.

    // -------------------------------------------------------------------------
    // Solver and linear algebra components
    // -------------------------------------------------------------------------
    mfem::CGSolver cgPP_solver; ///< Conjugate-gradient solver for φ_A.
    std::unique_ptr<mfem::HypreBoomerAMG> Mpp; ///< AMG preconditioner for stiffness matrix.

    // -------------------------------------------------------------------------
    // Linear system storage
    // -------------------------------------------------------------------------
    mfem::ParLinearForm B1t; ///< Linear form (reaction/forcing).
    mfem::HypreParVector X1v; ///< Solution vector (true DOFs).
    mfem::HypreParVector B1v; ///< RHS vector.
    mfem::HypreParVector Fpb; ///< Assembled forcing vector.
    mfem::HypreParVector Xs0; ///< Temporary solution vector.

    // -------------------------------------------------------------------------
    // Conductivity and reaction coefficients
    // -------------------------------------------------------------------------
    double gtPsi = 0.0; ///< Global total ψ (from Domain_Parameters).
    double gtPsA = 0.0; ///< Global anode ψ_A integral.

    mfem::ParGridFunction kap; ///< Conductivity κ_p.
    mfem::GridFunctionCoefficient cKp; ///< Coefficient wrapper for κ_p.
    mfem::GridFunctionCoefficient cRp; ///< Reaction coefficient.
    mfem::ParGridFunction RpP; ///< Reaction field in parallel space.
    mfem::ParGridFunction pP0; ///< Stored potential from previous iteration.

    // -------------------------------------------------------------------------
    // Bilinear and linear forms
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParBilinearForm> Kp2; ///< Conductivity-weighted stiffness form.
    std::unique_ptr<mfem::ParLinearForm>   Bp2; ///< Reaction linear form.
    mfem::HypreParMatrix KmP;                   ///< Stiffness matrix for φ_A.

    mfem::ParLinearForm Fpt; ///< Force term for particle potential.

    mfem::Array<int> ess_tdof_list_w; ///< List of essential true degrees of freedom for Dirichlet boundary conditions
    mfem::Array<int> dbc_w_bdr; ///< Array marking Dirichlet boundary attributes

};

#endif // POTA_HPP
