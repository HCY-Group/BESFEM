#pragma once

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "BoundaryConditions.hpp"
#include "Utils.hpp"
#include "FEMOperators.hpp"
#include <memory>

/**
 * @class PotentialBase
 * @brief Abstract base class for electrochemical potential solvers in BESFEM.
 *
 * This base class provides the common FEM infrastructure used for solving
 * electrochemical potentials in:
 * - **PotA** (anode potential)
 * - **PotC** (cathode potential)
 * - **PotE** (electrolyte potential)
 *
 * Responsibilities of PotentialBase:
 * - Store mesh/geometry and FE-space context
 * - Provide virtual interfaces for setup, assembly, and updates
 * - Manage temporary vectors and per-element volume info
 * - Serve as a uniform interface for potential solvers in the BESFEM pipeline
 *
 * Derived classes must implement:
 * - SetupField()       — Initialize potential field and assemble operators  
 * - AssembleSystem()   — Form system matrices/vectors based on C, ψ, etc.  
 * - UpdatePotential()  — Advance potential using reaction terms and FEM solves  
 */
class PotentialBase {
public:

    /**
     * @brief Construct the base potential solver.
     *
     * Stores geometry and domain parameters, initializes mesh pointers and
     * basic FEM quantities shared by PotA, PotC, and PotE.
     *
     * @param geo  Geometry/mesh handler (serial + parallel).
     * @param para Domain parameters (material constants, psi fields, etc.).
     */
    PotentialBase(Initialize_Geometry &geo, Domain_Parameters &para);

    /// Virtual destructor.
    virtual ~PotentialBase() = default;

    // -------------------------------------------------------------------------
    // Pure virtual interface implemented by derived classes
    // -------------------------------------------------------------------------

    /**
     * @brief Initialize the potential field and assemble required operators.
     *
     * Derived solvers may:
     * - initialize Dirichlet/Neumann boundary markers,
     * - set the starting potential,
     * - prepare preconditioners or mass/stiffness forms,
     * - mask invalid regions using phase-field ψ.
     *
     * @param ph            Potential field to initialize.
     * @param initial_value Scalar initial potential value.
     * @param psx           Phase-field mask ψ.
     */
    virtual void SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx) = 0;

    /**
     * @brief Assemble linear system operators for the current timestep.
     *
     * Typical responsibilities:
     * - build mass/stiffness contributions depending on concentration Cn,
     * - apply ψ-based masking,
     * - assemble RHS and boundary contributions,
     * - prepare Hypre matrices for solving.
     *
     * @param Cn        Concentration field used for conductivity terms.
     * @param psx       Phase-field mask ψ.
     * @param potential Potential field (input/output).
     */
    virtual void AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential) = 0;

    /**
     * @brief Solve for the updated potential and compute error metrics.
     *
     * Derived solvers compute:
     * - new potential phx via FEM/Hypre solve,
     * - reaction-related contributions,
     * - global L2/RMS error (returned in @p gerror).
     *
     * @param Rx     Reaction/source field.
     * @param phx    Potential field (input/output).
     * @param psx    Phase-field mask ψ.
     * @param gerror Output global error measure (MPI-reduced).
     */
    virtual void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) = 0;

    // -------------------------------------------------------------------------
    // General FEM state / diagnostics (inherited by all potential solvers)
    // -------------------------------------------------------------------------

    int nE = 0; ///< Number of elements.
    int nC = 0; ///< Nodes per element (corners).
    int nV = 0; ///< Total number of vertices.

    const mfem::Vector &EVol; ///< Per-element volumes from Domain_Parameters.

    Initialize_Geometry &geometry;          ///< Geometry and mesh context.
    Domain_Parameters   &domain_parameters; ///< Material/domain parameters.

    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh pointer.

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

    mfem::ParGridFunction TmpF; ///< Temporary field (e.g., residual/error workspace).

    mfem::ParGridFunction px0; ///< Current potential field (GridFunction).
    mfem::HypreParVector  X0;  ///< True DOF solution vector for the solver.
};

