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
 * potential, and reaction solvers. Particle pair arrays allocate only entries
 * [j][k] with j < k; reverse and diagonal entries remain null.
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
     * - Allocates phase masks, interface fields, and particle-pair storage
     * - Copies and clamps geometry masks, then constructs interfaces
     * - Computes element volumes (EVol)
     * - Integrates ψ and ψₑ to compute gtPsi, gtPse
     * - Computes global target current gTrgI
     * - Saves mesh and domain fields to the existing output_directory
     *
     */
    void SetupDomainParameters(const std::string& output_directory = ".");

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

    std::vector<double> gtPsA; ///< Global phase integrals in anode particle-group order.
    std::vector<double> gtPsC; ///< Global phase integrals in cathode particle-group order.

    mfem::Vector EVol; ///< Element volumes for FEM integration.

    std::unique_ptr<mfem::ParGridFunction> denom; ///< Denominator/workspace field used in phase-field normalization.

    std::unique_ptr<mfem::ParGridFunction> AvPA; ///< Magnitude of the total anode phase gradient.
    std::unique_ptr<mfem::ParGridFunction> AvPC; ///< Magnitude of the total cathode phase gradient.

    std::unique_ptr<mfem::ParGridFunction> denomA; ///< Sum of anode interface densities used to normalize weights.
    std::unique_ptr<mfem::ParGridFunction> denomC; ///< Sum of cathode interface densities used to normalize weights.

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

    std::vector<double> gTrgPsA; ///< Global per-particle anode target currents.
    std::vector<double> gTrgPsC; ///< Global per-particle cathode target currents.

    /// Reference to geometry handler.
    Initialize_Geometry &geometry;
    const SimulationConfig& cfg; ///< Borrowed runtime configuration; must outlive this object.


private:

    /// Owned fields for particle groups in configuration order.
    using FieldList = std::vector<std::unique_ptr<mfem::ParGridFunction>>;
    /// Pair fields indexed by [j][k]; only j < k entries are allocated.
    using PairFields = std::vector<FieldList>;

    /** @brief Borrowed view of one electrode's particle arrays.
     * Keeps the public half-cell/anode/cathode storage compatible with callers
     * while shared algorithms operate on the same set of named fields.
     */
    struct ParticleGroups
    {
        FieldList& phase; ///< Per-group phase masks.
        FieldList& gradient; ///< Per-group gradient magnitudes.
        FieldList& electrolyte_interface; ///< Electrolyte interface densities.
        FieldList& electrolyte_weight; ///< Electrolyte coupling weights.
        PairFields& pair_interface; ///< Particle-pair interface densities.
        PairFields& pair_phase; ///< Combined masks for each pair.
        PairFields& pair_weight; ///< Particle-pair coupling weights.
        std::vector<double>& local_total; ///< Local phase integrals.
        std::vector<double>& global_total; ///< MPI-reduced phase integrals.
        std::vector<double>& target_current; ///< MPI-reduced group target currents.
    };

    /** @brief Select the particle arrays for shared operations.
     * @param electrode ANODE or CATHODE; half cells use their single active set.
     * @return Non-owning references to this object's particle arrays.
     */
    ParticleGroups GetParticleGroups(sim::Electrode electrode);

    /** @brief Allocate group fields, upper-triangular pairs, and zero totals.
     * @param groups Borrowed particle arrays to populate.
     * @param count Number of particle groups.
     */
    void AllocateParticleGroups(ParticleGroups groups, std::size_t count);

    /** @brief Copy masks, accumulate their raw sum, and clamp individual masks.
     * @param groups Destination particle arrays.
     * @param masks Filtered geometry masks in group order.
     * @param[out] total Raw mask sum; the caller clamps it after combining electrodes.
     */
    void CopyParticleMasks(ParticleGroups groups, const FieldList& masks, mfem::ParGridFunction& total);

    /** @brief Build one electrode's interfaces, denominator, and coupling weights.
     * @param groups Particle masks and output interface storage.
     * @param[out] denominator Sum of pair and electrolyte interface densities.
     * @pre AvE contains the electrolyte gradient magnitude.
     */
    void BuildParticleInterfaces(ParticleGroups groups, mfem::ParGridFunction& denominator);

    /** @brief Integrate each group and calculate its material-dependent target current.
     * @param groups Phase masks and output totals.
     * @param materials One material per particle group.
     * @return Sum of the electrode's group target currents.
     * @pre EVol contains current mesh element volumes.
     */
    double CalculateParticleTotals(ParticleGroups groups, const std::vector<sim::MaterialType>& materials);


    // -------------------------------------------------------------------------
    // Internal setup routines
    // -------------------------------------------------------------------------

    /// Allocate phase and interface storage for the configured cell mode.
    void InitializeGridFunctions();

    /**
     * @brief Project/interpolate distance-function-based parameters.
     *
     * Copies the filtered masks from Initialize_Geometry, clamps phase values,
     * and builds interface fields and weights for the configured cell mode.
     *
     */
    void InterpolateDomainParameters();

    /**
     * @brief Allocate half-cell phase, interface, pair, and integration storage.
     */
    void InitializeHalfCellGridFunctions();
    /**
     * @brief Allocate separate anode and cathode phase, interface, and pair storage.
     */
    void InitializeFullCellGridFunctions();

    /**
     * @brief Copy and clamp filtered half-cell geometry masks into domain fields.
     */
    void InterpolateHalfCellMasks();
    /**
     * @brief Copy and clamp filtered electrode and electrolyte masks into domain fields.
     */
    void InterpolateFullCellMasks();

    /**
     * @brief Build half-cell gradient magnitudes, pair interfaces, and coupling weights.
     */
    void BuildHalfCellInterfaces();
    /**
     * @brief Build electrode-specific gradients, pair interfaces, and coupling weights.
     */
    void BuildFullCellInterfaces();

    /**
     * @brief Compute the Euclidean magnitude of a scalar field gradient.
     * @param phase_in Phase indicator.
     * @param[out] gradient_out Magnitude of its spatial gradient.
     */
    void ComputeGradientMagnitude(const mfem::ParGridFunction &phase_in, mfem::ParGridFunction &gradient_out);

    /** @brief Compute the geometric mean of two gradient-magnitude fields.
     * @param[out] out Interface density for a particle pair or electrolyte contact.
     * @param gradient_a First gradient magnitude.
     * @param gradient_b Second gradient magnitude.
     */
    void BuildInterface(mfem::ParGridFunction& out, const mfem::ParGridFunction& gradient_a,
        const mfem::ParGridFunction& gradient_b);

    /**
     * @brief Sum two phase masks and clamp their sum to [0, 1].
     * @param[out] out Combined pair mask.
     * @param phase_a First phase mask.
     * @param phase_b Second phase mask.
     */
    void BuildPairPhaseMask(mfem::ParGridFunction &out, const mfem::ParGridFunction &phase_a, const mfem::ParGridFunction &phase_b);

    /**
     * @brief Compute nonnegative interface-density ratios raised to power 0.8.
     * @param[out] weight_out Weights, zero where denominator is at most 1e-30.
     * @param numerator Interface density to normalize.
     * @param denominator Sum of interface densities.
     * @param mask Optional phase mask multiplied into the weights.
     */
    void ComputeInterfaceWeight(mfem::ParGridFunction &weight_out, const mfem::ParGridFunction &numerator, const mfem::ParGridFunction &denominator,
        const mfem::ParGridFunction *mask = nullptr);

    /**
     * @brief Integrate half-cell phase volumes and compute per-group target currents.
     */
    void CalculateHalfCellPhasePotentialsAndTargetCurrent();
    /**
     * @brief Integrate both electrodes and compute their group and total target currents.
     */
    void CalculateFullCellPhasePotentialsAndTargetCurrent();
    /**
     * @brief Convert a local phase volume into a globally reduced target current.
     * @param local_phase_volume Local integrated particle phase.
     * @param[out] global_target_current MPI-reduced target current.
     * @param material Material supplying lithium site density.
     */
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
     * @brief Integrate a phase field using cached element volumes.
     *
     * Calls CalculateTotals to compute the local and global integrals.
     * @pre EVol was filled by CalculatePhasePotentialsAndTargetCurrent().
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


    /**
     * @brief Print diagnostic totals (rank 0 only).
     *
     * Logs gtPsi, gtPse, gTrgI, and other totals for debugging or inspection.
     */
    void PrintInfo();

    // -------------------------------------------------------------------------
    // Geometry / storage members
    // -------------------------------------------------------------------------


    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh reference.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

    // -------------------------------------------------------------------------
    // Target values (set via integration)
    // -------------------------------------------------------------------------
    double tPsi = 0.0; ///< Local ψ total before MPI reduction.
    double tPse = 0.0; ///< Local ψₑ total before MPI reduction.

};

#endif // DOMAIN_PARAMETERS_HPP
