#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>
#include "SimTypes.hpp"

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class BoundaryConditions
 * @brief Handles all boundary-condition setup for BESFEM battery simulations.
 *
 * This class provides infrastructure for defining, marking, and imposing
 * Dirichlet/Neumann boundary conditions on the electrode and electrolyte
 * domains of a battery cell.
 *
 * The class stores both global and parallel mesh information, boundary
 * markers, and local-to-global element mappings required during FEM assembly.
 */
class BoundaryConditions {
public:

    /**
     * @brief Construct a BoundaryConditions object.
     *
     * Initializes geometry references, domain parameters, and internal data
     * structures required to construct the boundary masks. 
     *
     * @param geo  Reference to the initialized geometry (meshes, FE spaces, ψ-fields).
     * @param para Reference to global domain parameters (operating mode, material constants).
     */
    BoundaryConditions(Initialize_Geometry &geo, Domain_Parameters &para);

    /// Reference to simulation geometry data.
    Initialize_Geometry &geometry;

    /// Reference to domain parameters (operation mode, potentials, ψ-fields, etc.).
    Domain_Parameters &domain_parameters;

    /// Virtual destructor.
    virtual ~BoundaryConditions();

    /**
     * @brief Creates Dirichlet and Neumann boundary markers for the solver.
     *
     * This routine identifies West/East boundaries and electrolyte boundaries
     * based on geometry, mesh extents, and simulation mode (full cell vs half cell).
     *
     * @param mode      Cell mode (full-cell, anode-only, cathode-only).
     * @param electrode Active electrode domain (anode or cathode).
     */
    void SetupBoundaryConditions(sim::CellMode mode, sim::Electrode electrode);

    /**
     * @brief Sets the pinned DOF for electrolyte potential stability.
     *
     * Solving for the electrolyte potential for the full cell requires anchoring a single DOF
     * to avoid a null space. This routine selects a vertex inside the
     * electrolyte ψ-field (above threshold) and far from domain boundaries.
     *
     */
    void SetupPinnedDOF(mfem::ParFiniteElementSpace &fespace);

    /**
     * @brief Selects a pinned node closest to the geometric center of the electrolyte.
     *
     * Searches all elements owned by the current MPI rank, selects a vertex
     * whose ψ_E value exceeds the threshold and whose distance to the
     * domain center is minimal.
     *
     * @param threshold Minimum ψ_E value to consider a point inside the electrolyte.
     * @param out_dist2 Output: squared distance of the selected point to the center.
     * @return Global vertex index of the pinned node, or -1 if none found.
     */
    int SelectCenterPin(double threshold, double &out_dist2);

    /**
     * @brief Returns the first acceptable pinned node found on the current rank.
     *
     * Similar to SelectCenterPin(), but stops early at the first vertex that:
     * - is inside the electrolyte (ψ_E ≥ threshold),
     * - is not too close to West/East boundaries,
     * - lies away from min/max Z.
     *
     * Useful when fast selection is needed.
     *
     * @param threshold Minimum ψ_E value to consider a vertex inside the electrolyte.
     * @param out_dist2 Output: squared distance of chosen vertex to center.
     * @return Global vertex index of the selected DOF, or -1 if none found.
     */
    int SelectFirstPin(double threshold, double &out_dist2);


    // -------------------------------------------------------------------------
    // PUBLIC MEMBERS (Boundary markers, DOF lists, mesh storage)
    // -------------------------------------------------------------------------

    int gei; ///< Global element index (used during element iteration).

    mfem::Array<int> nbc_w_bdr; ///< Neumann boundary markers for West boundary.
    mfem::Array<int> nbc_e_bdr; ///< Neumann boundary markers for East boundary.
    mfem::Array<int> nbc_bdr;   ///< Neumann markers for electrolyte boundaries.

    mfem::Array<int> dbc_bdr;   ///< Dirichlet markers for electrolyte boundaries.
    mfem::Array<int> dbc_w_bdr; ///< Dirichlet markers for West boundary.
    mfem::Array<int> dbc_e_bdr; ///< Dirichlet markers for East boundary.

    mfem::Array<int> ess_tdof_list_w; ///< Essential DOFs on West boundary.
    mfem::Array<int> ess_tdof_list_e; ///< Essential DOFs on East boundary.

    mfem::Array<int> gVTX; ///< Global vertex indices of important geometry vertices (corners).
    mfem::Array<int> VTX;  ///< Local vertex indices of corners.

    /// Global serial mesh (rank 0 only).
    mfem::Mesh globalMesh;

    /// Parallel distributed mesh.
    mfem::ParMesh parallelMesh;

    /// Finite element space on the parallel mesh.
    mfem::ParFiniteElementSpace parfespace;

    /// Local-to-global element index mapping.
    mfem::Array<HYPRE_BigInt> E_L2G;

    bool pin; ///< True if a pinned DOF has been assigned for the electrolyte potential.

    mfem::Array<int> ess_tdof_list; ///< Full list of essential DOFs for Dirichlet BCs.
    mfem::Array<int> ess_tdof_marker; ///< Marker array for essential (Dirichlet) boundaries.

    int myid; ///< MPI rank ID.
    int rkpp; ///< MPI rank that owns the pinned DOF.
};

#endif // BOUNDARYCONDITIONS_HPP
