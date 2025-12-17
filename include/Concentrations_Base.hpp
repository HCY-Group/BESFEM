#pragma once

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "BoundaryConditions.hpp"
#include "Utils.hpp"
#include "FEMOperators.hpp"
#include <memory>

/**
 * @class ConcentrationBase
 * @brief Abstract base class for all concentration solvers in BESFEM.
 *
 * Provides the common interface and shared data structures used by all
 * concentration-update modules:
 * - **CnA**  (Cahn–Hilliard solid anode)
 * - **CnC**  (diffusion-based cathode)
 * - **CnE**  (electrolyte diffusion + reaction coupling)
 *
 * Derived classes must implement:
 * - SetupField()  → initialize concentration & assemble matrices  
 * - UpdateConcentration() → advance concentration per timestep  
 *
 * ## Responsibilities of ConcentrationBase
 * - Store mesh, geometry, element volumes, and FE-space references  
 * - Track global lithiation/salt fractions (Xfr, CeC, gCeC, CeAvg)  
 * - Provide workspace buffers for per-element averaging, masks, and scratch vectors  
 * - Serve as a shared foundation for MFEM/Hypre-based concentration solvers  
 */
class ConcentrationBase {
public:

    /**
     * @brief Construct the base concentration handler.
     *
     * Stores references to geometry and domain parameters, allocates shared
     * buffers (element averages, vertex buffer), and initializes FE data.
     *
     * @param geo  Geometry handler (mesh, FE space).
     * @param para Domain parameters (material constants, element volumes, etc.).
     */
    ConcentrationBase(Initialize_Geometry &geo, Domain_Parameters &para);

    /// Virtual destructor.
    virtual ~ConcentrationBase() = default;

    /**
     * @brief Initialize the concentration field and assemble operator components.
     *
     * Derived classes must implement model-specific assembly:
     * - mass/stiffness matrices
     * - reaction/mobility/diffusivity fields
     * - initial condition masking
     *
     * @param Cn            Concentration field to initialize.
     * @param initial_value Scalar initial value for `Cn`.
     * @param psx           Phase-field ψ used for region masking.
     */
    virtual void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx) = 0;

    /**
     * @brief Advance the concentration field by one timestep.
     *
     * Each derived solver implements its own update method:
     * - CnA → Cahn–Hilliard  
     * - CnC → standard diffusion  
     * - CnE → electrolyte CN update  
     *
     * @param Rx  Reaction/source field for this timestep.
     * @param Cn  Concentration field (input/output).
     * @param psx Phase-field ψ for masking region.
     */
    virtual void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn,
                                     mfem::ParGridFunction &psx) = 0;

    /**
     * @brief Return the current degree of lithiation or salt fraction.
     *
     * The value depends on the solver:
     * - CnA → anode lithium fraction  
     * - CnC → cathode lithium fraction  
     * - CnE → electrolyte salt inventory  
     */
    double GetLithiation() const { return Xfr; }

    // -------------------------------------------------------------------------
    // Global state / diagnostics
    // -------------------------------------------------------------------------
    double Xfr  = 0.0; ///< Global degree of lithiation or concentration.
    double CeC  = 0.0; ///< Local (per-element) salt or concentration moment.
    double gCeC = 0.0; ///< Global salt inventory (summed over domain).

    double CeAvg = 0.0; ///< Global average concentration.
    double Ce0   = 0.001; ///< Initial background concentration (reference).

    // -------------------------------------------------------------------------
    // Mesh / element / vertex information
    // -------------------------------------------------------------------------
    mfem::Array<double> VtxVal; ///< Vertex-value workspace buffer.

    int nE = 0; ///< Number of mesh elements.
    int nC = 0; ///< Number of nodes per element (corners).
    int nV = 0; ///< Total number of global vertices.

    mfem::Vector EAvg;          ///< Per-element averages (workspace).
    const mfem::Vector &EVol;   ///< Element volumes (from Domain_Parameters).

    // -------------------------------------------------------------------------
    // Phase-field normalization
    // -------------------------------------------------------------------------
    double gtPsi = 0.0; ///< Global ψ normalization (total).
    double gtPse = 0.0; ///< Global ψ normalization (electrolyte or solid).

    // -------------------------------------------------------------------------
    // Boundary and workspace DoF storage
    // -------------------------------------------------------------------------
    mfem::Array<int> boundary_dofs; ///< Boundary true DOFs.
    mfem::HypreParVector X1v;       ///< Scratch vector for assembly/solves.

    // -------------------------------------------------------------------------
    // Geometry and FE context
    // -------------------------------------------------------------------------
    Initialize_Geometry &geometry;          ///< Geometry and mesh context.
    Domain_Parameters   &domain_parameters; ///< Physics/material data.

    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh (distributed).
    mfem::Mesh    *gmesh = nullptr; ///< Global serial mesh (rank 0, if available).

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

    std::unique_ptr<mfem::ParGridFunction> TmpF; ///< Temporary grid-function buffer.
};

