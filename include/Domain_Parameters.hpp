#ifndef DOMAIN_PARAMETERS_HPP
#define DOMAIN_PARAMETERS_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>

/**
 * @file Domain_Parameters.hpp
 * @brief Defines the Domain_Parameters class for initializing and managing
 *        domain-specific fields and integrals in BESFEM battery simulations.
 *
 * This class constructs and stores phase fields (ψ, ψₑ, ψ_A, ψ_C),
 * surface-area fields (AvP, AvA, AvC, AvB), distance functions, and
 * element-volume data. It also computes global integrals such as gtPsi,
 * gtPse, gtPsA, gtPsC, and the global target current gTrgI.
 */

class Initialize_Geometry;

/**
 * @class Domain_Parameters
 * @brief Manages domain-specific phase fields, auxiliary fields, and global integrals.
 *
 * The Domain_Parameters class provides:
 * - Construction and interpolation of phase fields: ψ (solid), ψₑ (electrolyte), ψ_A, ψ_C  
 * - Computation of auxiliary “surface-area density” fields AvP, AvA, AvC, AvB  
 * - Element volumes (EVol) for FEM integrations  
 * - Computation of global totals (gtPsi, gtPse, gtPsA, gtPsC)  
 * - Calculation of the global target current (gTrgI) used for current-controlled simulations  
 *
 * The results are used by concentration solvers (CnA, CnC, CnE) and potential solvers.
 */
class Domain_Parameters {

public:

    /**
     * @brief Construct a Domain_Parameters object.
     *
     * Stores a reference to geometry, allocates grid functions,
     * initializes internal counters (nE, nV, nC), and prepares
     * distance-function pointers.
     *
     * @param geo Reference to the initialized geometry (mesh, FE space, ψ distance fields).
     */
    Domain_Parameters(Initialize_Geometry &geo);

    /// Destructor.
    virtual ~Domain_Parameters();

    /**
     * @brief Initialize all domain parameters based on the mesh type.
     *
     * Steps:
     * - Allocates and zeros ψ, ψₑ, ψ_A, ψ_C, AvP, AvA, AvC, AvB  
     * - Projects/interpolates distance-function-based fields  
     * - Computes element volumes (EVol)  
     * - Integrates ψ and ψₑ to compute gtPsi, gtPse  
     * - Computes global target current gTrgI  
     *
     * @param mesh_type Character flag identifying the mesh geometry:
     *        - `"ml"` MATLAB mesh  
     *        - `"v"` voxel-derived mesh  
     */
    void SetupDomainParameters(const char* mesh_type);

    // -------------------------------------------------------------------------
    // Phase fields (grid functions)
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParGridFunction> psi; ///< Solid-phase indicator (ψ).
    std::unique_ptr<mfem::ParGridFunction> pse; ///< Electrolyte-phase indicator (ψₑ).
    std::unique_ptr<mfem::ParGridFunction> psA; ///< Anode-phase indicator.
    std::unique_ptr<mfem::ParGridFunction> psC; ///< Cathode-phase indicator.

    std::unique_ptr<mfem::ParGridFunction> denom;

    // -------------------------------------------------------------------------
    // Surface-area / geometry-related auxiliary fields
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParGridFunction> AvP; ///< Particle surface-area density.
    std::unique_ptr<mfem::ParGridFunction> AvA; ///< Anode surface-area density.
    std::unique_ptr<mfem::ParGridFunction> AvC; ///< Cathode surface-area density.
    std::unique_ptr<mfem::ParGridFunction> AvB; ///< Boundary surface-area density.
    std::unique_ptr<mfem::ParGridFunction> AvE; ///< Electrolyte surface-area density.

    // -------------------------------------------------------------------------
    // Global integrals and target current
    // -------------------------------------------------------------------------
    double gtPsi = 0.0; ///< Global integral of ψ (solid).
    double gtPse = 0.0; ///< Global integral of ψₑ (electrolyte).
    double gTrgI = 0.0; ///< Global target current (galvanostatic control).

    double gTrg1 = 0.0; ///< Global target current for phase 1.
    double gTrg2 = 0.0; ///< Global target current for phase 2.
    double gTrg3 = 0.0; ///< Global target current for phase 3.

    double gtPsi1 = 0.0; ///< Global integral of ψ for phase 1.
    double gtPsi2 = 0.0; ///< Global integral of ψ for phase
    double gtPsi3 = 0.0; ///< Global integral of ψ for phase 3.

    double gtPsA = 0.0; ///< Global integral of ψ_A (anode region).
    double gtPsC = 0.0; ///< Global integral of ψ_C (cathode region).

    mfem::Vector EVol; ///< Element volumes for FEM integration.

    std::vector<int> particle_labels;
    std::vector<std::unique_ptr<mfem::ParGridFunction>> ps;
    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvPs;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> AvP_Pairs;
    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvEs;
    std::vector<std::unique_ptr<mfem::ParGridFunction>> WeightEs;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> psi_Pairs;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> WeightPairs;
    std::vector<double> tPs;
    std::vector<double> gtPs;
    std::vector<double> gTrgPs;

    /// Reference to geometry handler.
    Initialize_Geometry &geometry;

private:

    // -------------------------------------------------------------------------
    // Internal setup routines
    // -------------------------------------------------------------------------

    /// Allocate grid functions (psi, pse, AvP, AvA, AvC, AvB).
    void InitializeGridFunctions();

    /**
     * @brief Project/interpolate distance-function-based parameters.
     *
     * Converts dsF / dsF_A / dsF_C into phase indicators ψ, ψₑ and surface-area
     * fields AvP, etc., depending on mesh type. Clamps values to ensure
     * physical consistency.
     *
     * @param mesh_type Mesh type specifier ("ml", "v").
     */
    void InterpolateDomainParameters(const char* mesh_type);

    /**
     * @brief Compute the local and global totals of a field.
     *
     * Performs:
     * - Element-wise multiplication with the element volumes (EVol)  
     * - Local summation  
     * - MPI reduction to compute a global total  
     *
     * @param grid_function Field to integrate.
     * @param element_volumes Precomputed element volumes.
     * @param local_total [out] Process-local integral.
     * @param global_total [out] MPI-reduced total.
     */
    void CalculateTotals(const mfem::ParGridFunction &grid_function, const mfem::Vector &element_volumes,
                         double &local_total, double &global_total);

    /**
     * @brief Compute EVol and integrate a phase field (ψ or ψₑ).
     *
     * Fills @ref EVol using geometric data, then calls @ref CalculateTotals
     * to compute the local and global integrals.
     *
     * @param grid_function Phase-field indicator.
     * @param total         [out] Local total.
     * @param global_total  [out] Global total (MPI).
     */
    void CalculateTotalPhaseField(const mfem::ParGridFunction &grid_function, double &total,
                                  double &global_total);

    /**
     * @brief Compute ψ, ψₑ integrals and derive the global target current.
     *
     * Calls @ref CalculateTotalPhaseField for ψ and ψₑ, then calls
     * @ref CalculateTargetCurrent to update gTrgI.
     */
    void CalculatePhasePotentialsAndTargetCurrent();

    /**
     * @brief Compute local contribution to target current from ψ.
     *
     * Uses electrochemical constants (density, capacity, etc.) to convert
     * the integrated particle-phase volume into a target current contribution.
     *
     * @param total_psi Local integral of ψ.
     * @param global_total [out] Global target current contribution after MPI reduction.
     */
    void CalculateTargetCurrent(double total_psi, double &global_total);

    /**
     * @brief Print diagnostic totals (rank 0 only).
     *
     * Logs gtPsi, gtPse, gTrgI, and other totals for debugging or inspection.
     */
    void PrintInfo();

    // -------------------------------------------------------------------------
    // Geometry / storage members
    // -------------------------------------------------------------------------
    int nV = 0; ///< Number of vertices.
    int nE = 0; ///< Number of elements.
    int nC = 0; ///< Nodes per element (corners).

    mfem::ParGridFunction *dsF     = nullptr; ///< Distance function (whole domain).
    mfem::ParGridFunction *dsF_A   = nullptr; ///< Distance function (anode).
    mfem::ParGridFunction *dsF_C   = nullptr; ///< Distance function (cathode).

    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh reference.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

    // -------------------------------------------------------------------------
    // Target values (set via integration)
    // -------------------------------------------------------------------------
    double tPsi = 0.0; ///< Local ψ total before MPI reduction.
    double tPse = 0.0; ///< Local ψₑ total before MPI reduction.
    double trgI = 0.0; ///< Local target current before global reduction.

    double tPsi1 = 0.0; ///< Local ψ for phase 1 total before MPI reduction.
    double tPsi2 = 0.0; ///< Local ψ for phase 2 total before MPI reduction.
    double tPsi3 = 0.0; ///< Local ψ for phase 3 total before MPI reduction.

    double tPsA = 0.0; ///< Local ψ_A total before MPI reduction.
    double tPsC = 0.0; ///< Local ψ_C total before MPI reduction.
};

#endif // DOMAIN_PARAMETERS_HPP
