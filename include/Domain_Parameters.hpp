#ifndef DOMAIN_PARAMETERS_HPP
#define DOMAIN_PARAMETERS_HPP

#include "mfem.hpp"
#include "SimulationConfig.hpp"
#include <memory>
#include <string>
#include <vector>

/**
 * @file Domain_Parameters.hpp
 * @brief Defines domain fields, geometry fields, and global integrals used by BESFEM.
 */

class Initialize_Geometry;

/**
 * @class Domain_Parameters
 * @brief Stores geometry-dependent fields and global quantities used throughout BESFEM.
 *
 * Domain_Parameters constructs and manages the phase-field masks, interface
 * fields, element volumes, and global integrals required by the concentration,
 * potential, and reaction solvers.
 */
class Domain_Parameters {

public:

    /**
     * @brief Construct a Domain_Parameters object.
     *
     * Stores references to the geometry and simulation configuration, initializes
     * mesh/finite-element-space pointers, and prepares storage for domain fields.
     *
     * @param geo Reference to the initialized geometry object.
     * @param cfg Reference to the simulation configuration.
     */
    Domain_Parameters(Initialize_Geometry &geo, const SimulationConfig &cfg);

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
     */
    void SetupDomainParameters();

    // -------------------------------------------------------------------------
    // Phase fields (grid functions)
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParGridFunction> psi; ///< Solid-phase indicator (ψ).
    std::unique_ptr<mfem::ParGridFunction> pse; ///< Electrolyte-phase indicator (ψₑ).
    std::unique_ptr<mfem::ParGridFunction> psiA; ///< Anode-phase indicator.
    std::unique_ptr<mfem::ParGridFunction> psiC; ///< Cathode-phase indicator.

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
   
    double gtPsiA = 0.0; ///< Global integral of ψ_A (anode).
    double gtPsiC = 0.0; ///< Global integral of ψ_C (cathode).

    double tPsiA = 0.0; ///< Local integral of ψ_A (anode).
    double tPsiC = 0.0; ///< Local integral of ψ_C (cathode).   

    double gTrgIA = 0.0; ///< Global target current for anode.
    double gTrgIC = 0.0; ///< Global target current for cathode.

    double gTrg1 = 0.0; ///< Global target current for phase 1.
    double gTrg2 = 0.0; ///< Global target current for phase 2.
    double gTrg3 = 0.0; ///< Global target current for phase 3.

    double gtPsi1 = 0.0; ///< Global integral of ψ for phase 1.
    double gtPsi2 = 0.0; ///< Global integral of ψ for phase 2.
    double gtPsi3 = 0.0; ///< Global integral of ψ for phase 3.


    std::vector<double> gtPsA;
    std::vector<double> gtPsC;

    mfem::Vector EVol; ///< Element volumes for FEM integration.

    std::unique_ptr<mfem::ParGridFunction> denom; ///< Denominator/workspace field used in phase-field normalization.

    std::unique_ptr<mfem::ParGridFunction> AvPA;
    std::unique_ptr<mfem::ParGridFunction> AvPC;

    std::unique_ptr<mfem::ParGridFunction> denomA;
    std::unique_ptr<mfem::ParGridFunction> denomC;

    std::vector<int> particle_labels; ///< Material/particle labels read from the segmented geometry.
    std::vector<int> anode_particle_labels; ///< Anode particle labels.
    std::vector<int> cathode_particle_labels; ///< Cathode particle labels.

    std::vector<std::unique_ptr<mfem::ParGridFunction>> ps; ///< Per-particle phase-field masks.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvPs; ///< Per-particle surface-area density fields.

    std::vector<std::unique_ptr<mfem::ParGridFunction>> psA; ///< Per-particle phase fields anode.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> psC; ///< Per-particle phase fields cathode.

    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvPsA; ///< Per-particle surface-area density fields anode.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvPsC; ///< Per-particle surface-area density fields cathode.

    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> AvP_Pairs; ///< Pairwise particle-particle interfacial area fields.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> AvP_PairsA; ///< Pairwise anode particle-particle interfacial area fields.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> AvP_PairsC; ///< Pairwise cathode particle-particle interfacial area fields.   

    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvEs; ///< Per-particle electrolyte interface area fields.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvEsA; ///< Per-particle anode electrolyte interface area fields.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> AvEsC; ///< Per-particle cathode electrolyte interface area fields.

    std::vector<std::unique_ptr<mfem::ParGridFunction>> WeightEs; ///< Per-particle electrolyte coupling weights.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> WeightEsA; ///< Per-particle anode electrolyte coupling weights.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> WeightEsC; ///< Per-particle cathode electrolyte coupling weights.

    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> psi_Pairs; ///< Pairwise particle-particle interface phase fields.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> psi_PairsA; ///< Pairwise anode particle-particle interface phase fields.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> psi_PairsC; ///< Pairwise cathode particle-particle interface phase fields.

    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> WeightPairs; ///< Pairwise particle-particle coupling weights.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> WeightPairsA; ///< Pairwise anode particle-particle coupling weights.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> WeightPairsC; ///< Pairwise cathode particle-particle coupling weights.

    std::vector<double> tPs; ///< Local per-particle phase-field totals.
    std::vector<double> gtPs; ///< Global per-particle phase-field totals.
    std::vector<double> gTrgPs; ///< Global per-particle target currents.

    std::vector<double> tPsA; ///< Local per-particle anode phase-field totals.
    std::vector<double> tPsC; ///< Local per-particle cathode phase-field totals.

    // std::vector<double> gtPsA; ///< Global per-particle anode phase-field totals.
    // std::vector<double> gtPsC; ///< Global per-particle cathode phase-field totals.

    std::vector<double> gTrgPsA; ///< Global per-particle anode target currents.
    std::vector<double> gTrgPsC; ///< Global per-particle cathode target currents.

    /// Reference to geometry handler.
    Initialize_Geometry &geometry;
    const SimulationConfig& cfg;


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
     */
    void InterpolateDomainParameters();

    void InitializeHalfCellGridFunctions();
    void InitializeFullCellGridFunctions();

    void InterpolateHalfCellMasks();
    void InterpolateFullCellMasks();

    void BuildHalfCellInterfaces();
    void BuildFullCellInterfaces();

    // void ApplyAMR();

    void ComputeGradientMagnitude(const mfem::ParGridFunction &phase_in, mfem::ParGridFunction &gradient_out);

    void BuildPairInterface(mfem::ParGridFunction &out, const mfem::ParGridFunction &phase_a, const mfem::ParGridFunction &phase_b,
        const mfem::ParGridFunction &gradient_a, const mfem::ParGridFunction &gradient_b);

    void BuildElectrolyteInterface(mfem::ParGridFunction &out, const mfem::ParGridFunction &electrolyte_phase, const mfem::ParGridFunction &particle_gradient);
    void BuildPairPhaseMask(mfem::ParGridFunction &out, const mfem::ParGridFunction &phase_a, const mfem::ParGridFunction &phase_b);

    void ComputeInterfaceWeight(mfem::ParGridFunction &weight_out, const mfem::ParGridFunction &numerator, const mfem::ParGridFunction &denominator,
        const mfem::ParGridFunction *mask = nullptr);


    void CalculateHalfCellPhasePotentialsAndTargetCurrent();
    void CalculateFullCellPhasePotentialsAndTargetCurrent();
    void CalculateTargetCurrent(double local_phase_volume, double &global_target_current, sim::MaterialType material);

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
    void CalculateTotals(const mfem::ParGridFunction &grid_function, const mfem::Vector &element_volumes, double &local_total, double &global_total);

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
    void CalculateTotalPhaseField(const mfem::ParGridFunction &grid_function, double &total, double &global_total);

    /**
     * @brief Compute phase-field integrals and target currents.
     *
     * Computes total phase-field weights for solid, electrolyte, anode, cathode,
     * and per-particle fields, then updates the corresponding target-current
     * values.
     */
    void CalculatePhasePotentialsAndTargetCurrent();

    // /**
    //  * @brief Compute local contribution to target current from ψ.
    //  *
    //  * Uses electrochemical constants (density, capacity, etc.) to convert
    //  * the integrated particle-phase volume into a target current contribution.
    //  *
    //  * @param total_psi Local integral of ψ.
    //  * @param global_total [out] Global target current contribution after MPI reduction.
    //  */
    // void CalculateTargetCurrent(double total_psi, double &global_total, sim::MaterialType material);

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

    // double tPsA = 0.0; ///< Local ψ_A total before MPI reduction.
    // double tPsC = 0.0; ///< Local ψ_C total before MPI reduction.
};

#endif // DOMAIN_PARAMETERS_HPP
