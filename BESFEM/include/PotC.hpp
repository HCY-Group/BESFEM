#ifndef POTC_HPP
#define POTC_HPP

/**
 * @file PotC.hpp
 * @brief Cathode potential solver for BESFEM.
 *
 * Implements initialization, conductivity assembly, and time-stepping for
 * the solid active-material (particle) potential field φ_C. This class
 * extends the generic PotentialBase interface to handle cathode-specific
 * boundary conditions, conductivity fields, and reaction coupling.
 */

#include "Potentials_Base.hpp"

/**
 * @brief Boundary value for the cathode solid-phase potential φ_C.
 */
extern double BvC;

/**
 * @class PotC
 * @brief Solid-phase cathode potential solver.
 *
 * Responsibilities:
 * - Initialize φ_C and apply cathode boundary conditions
 * - Assemble conductivity-weighted stiffness operators
 * - Assemble reaction/forcing terms
 * - Solve the linear system each timestep (Hypre CG + AMG)
 * - Compute global error metrics
 */
class PotC : public PotentialBase {

public:

    /**
     * @brief Construct the cathode potential solver.
     *
     * Stores references to geometry, domain parameters, and boundary
     * conditions and creates helper classes for FEM assembly.
     *
     * @param geo  Geometry/mesh handler (serial + parallel).
     * @param para Domain parameter container.
     * @param bc   Boundary condition handler.
     */
    PotC(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc);

    Initialize_Geometry  &geometry;            ///< Geometry and mesh infrastructure
    Domain_Parameters    &domain_parameters;   ///< Material/phase-field parameters
    BoundaryConditions   &boundary_conditions; ///< Electrode boundary configuration

    FEMOperators fem; ///< FEM operator assembly utilities
    Utils utils;      ///< Utility helpers (clamping, reductions, etc.)

    double BvC;       ///< Cathode Dirichlet boundary value for φ_C (if used)

    /**
     * @brief Initialize solid-phase potential φ_C.
     *
     * Sets the initial potential, configures essential DOFs, assembles
     * conductivity operators, and initializes solver/preconditioner.
     *
     * @param ph            Potential field to initialize.
     * @param initial_value Scalar initial potential value.
     * @param psx           Phase-field mask ψ_C.
     */
    void SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx) override;

    /**
     * @brief Assemble the linear system for φ_C.
     *
     * Builds:
     * - conductivity-weighted stiffness matrix,
     * - reaction linear form,
     * - boundary condition modifications.
     *
     * @param Cn        Solid concentration field.
     * @param psx       Phase-field mask ψ_C.
     * @param potential Potential φ_C grid function (in/out).
     */
    void AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                        mfem::ParGridFunction &potential) override;

    /**
     * @brief Advance φ_C by one timestep.
     *
     * Recomputes reaction terms, updates the forcing vector, and solves the
     * cathode potential system. Computes a global error metric.
     *
     * @param Rx     Reaction source term.
     * @param phx    Potential field φ_C (in/out).
     * @param psx    Phase-field mask ψ_C.
     * @param gerror Output: global MPI-reduced RMS/L2 error metric.
     */
    void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx,
                         mfem::ParGridFunction &psx, double &gerror) override;

private:

    /**
     * @brief Compute solid-phase conductivity κ_p(C, ψ_C).
     *
     * Populates the conductivity field based on material rules and ψ_C masking.
     *
     * @param Cn  Concentration field.
     * @param psx Phase-field mask ψ_C.
     */
    void ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space for φ_C

    // -------------------------------------------------------------------------
    // Solver and linear algebra components
    // -------------------------------------------------------------------------
    mfem::CGSolver cgPP_solver;                     ///< Conjugate-gradient solver
    std::unique_ptr<mfem::HypreBoomerAMG> Mpp;      ///< AMG preconditioner

    // -------------------------------------------------------------------------
    // Linear system storage
    // -------------------------------------------------------------------------
    mfem::ParLinearForm B1t;     ///< Linear form for reaction/forcing
    mfem::HypreParVector X1v;    ///< Solution vector (true DOFs)
    mfem::HypreParVector B1v;    ///< RHS vector
    mfem::HypreParVector Fpb;    ///< Assembled forcing vector
    mfem::HypreParVector Xs0;    ///< Temporary solution vector

    // -------------------------------------------------------------------------
    // Conductivity and reaction coefficients
    // -------------------------------------------------------------------------
    double gtPsi = 0.0; ///< Global ψ from Domain_Parameters
    double gtPsC = 0.0; ///< Global ψ_C integral

    mfem::ParGridFunction kap;   ///< Conductivity κ_p(C, ψ_C)
    mfem::GridFunctionCoefficient cKp; ///< Conductivity coefficient wrapper
    mfem::GridFunctionCoefficient cRp; ///< Reaction coefficient wrapper
    mfem::ParGridFunction RpP;   ///< Reaction mapped into parallel space
    mfem::ParGridFunction pP0;   ///< Stored previous potential

    // -------------------------------------------------------------------------
    // Bilinear and linear forms
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParBilinearForm> Kp2; ///< Conductivity-weighted stiffness form
    std::unique_ptr<mfem::ParLinearForm>   Bp2; ///< Reaction linear form
    mfem::HypreParMatrix KmP;                   ///< Stiffness matrix for φ_C

    mfem::ParLinearForm Fpt; ///< Force term for particle potential

    mfem::Array<int> ess_tdof_list_e; ///< List of essential true degrees of freedom for Dirichlet boundary conditions
    mfem::Array<int> dbc_e_bdr; ///< Array marking Dirichlet boundary attributes


};

#endif // POTC_HPP
