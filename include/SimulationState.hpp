/**
 * @file SimulationState.hpp
 * @brief Runtime fields, solver ownership, and electrode update helpers.
 *
 * SimulationState owns the simulation data; battery_simulation.cpp controls
 * the timestep loop and calls the initialization and update methods below.
 */
#ifndef SIMULATION_STATE_HPP
#define SIMULATION_STATE_HPP

#include "mfem.hpp"
#include "BESFEM_All.hpp"

#include <memory>
#include <vector>

/**
 * @brief Concentration and reaction state for one particle group.
 *
 * Each group retains its own material, fields, and solvers, even when several
 * groups use the same material. Its electrode supplies the shared potential.
 */
struct ParticleState
{
    int label = -1; ///< Mesh/domain label identifying this particle group.
    sim::MaterialType material = sim::MaterialType::Graphite; ///< Material used for this group's transport and reaction properties.
    std::unique_ptr<ConcentrationBase> concentration; ///< Material-specific concentration solver.
    std::unique_ptr<mfem::ParGridFunction> Cn_gf; ///< Particle-group concentration field.
    std::unique_ptr<mfem::ParGridFunction> Cn_gf_psi; ///< Workspace for concentration masked by the particle phase field.
    std::unique_ptr<Reaction> reaction; ///< Owned reaction helper.
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf; ///< Reaction field; per-group in ParticleState, combined in half-cell SimulationState.
    std::unique_ptr<mfem::ParGridFunction> Rx_src; ///< Frozen reaction source used during a concentration update.

    /**
     * @brief Allocate and initialize one group's concentration and reaction state.
     *
     * Graphite and LFP use ElectrodeCahnHilliard; Carbon, Silicon, and NMC
     * use ElectrodeDiffusion. The reaction field starts at Constants::init_Rxn.
     * @param geometry Mesh and finite-element infrastructure.
     * @param domain_parameters Domain masks and global phase totals.
     * @param cfg Simulation configuration.
     * @param particle_material Material assigned to this group.
     * @param particle_label Domain label identifying the group.
     * @param init_cn Initial concentration for this group.
     * @param particle_field Spatial phase mask for this group.
     * @param particle_total Global phase total for this group.
     */
    void Initialize(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        const SimulationConfig& cfg, sim::MaterialType particle_material, int particle_label,
        double init_cn, mfem::ParGridFunction& particle_field, double particle_total);
};

/**
 * @brief Chemical-potential and concentration workspaces for particle pairs.
 *
 * Only entries [j][k] with j < k are allocated. A pair is stored once;
 * BuildPairTerms() selects the correct self and neighbor orientation.
 */
struct PairWorkspaces
{
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_a; ///< Chemical potential of the lower-index group in each pair.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_b; ///< Chemical potential of the higher-index group in each pair.
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> sum_pairs; ///< Pair concentration-sum workspaces, indexed by [j][k] with j < k.

    /**
     * @brief Allocate and zero workspaces for all np * (np - 1) / 2 pairs.
     * @param geometry Finite-element infrastructure used to allocate fields.
     * @param np Number of particle groups in this electrode.
     * @param electrode_name Electrode name used in diagnostic output.
     */
    void Initialize(Initialize_Geometry& geometry, int np, const char* electrode_name);

    /**
     * @brief Rebuild the pair-coupling list for particle group j.
     *
     * Returned entries borrow pointers to the pair fields; ownership stays with
     * these workspaces and the supplied domain fields.
     * @param weight_pairs Pair coupling weights, indexed by the ordered pair.
     * @param avp_pairs Pair interface fields, indexed by the ordered pair.
     * @param j Index of the group whose concentration will be updated.
     * @param[out] pair_terms Cleared and filled with one entry per other group.
     * @param np Number of particle groups.
     * @pre Required pair workspaces and domain fields are allocated.
     */
    void BuildPairTerms(const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs, int j,
        std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np) const;
};

/**
 * @brief State and update helpers for one electrode.
 *
 * Owns one spatially varying solid potential field shared by all particle
 * groups. Anode and cathode instances own separate fields and solvers.
 */
struct ElectrodeState
{
    std::vector<ParticleState> particles; ///< One state per particle group, in configuration/domain group order.
    PairWorkspaces pairs; ///< Workspaces for coupling groups within this electrode.
    std::unique_ptr<ElectrodePotential> potential; ///< Shared solid-phase potential solver for this electrode.
    std::unique_ptr<mfem::ParGridFunction> ph_gf; ///< Spatially varying solid-phase potential field.

    /**
     * @brief Initialize the shared potential, particle states, and pair workspaces.
     *
     * Uses the first material for initial potential setup. Subsequent multi-group
     * potential assembly receives all materials through BuildParticleFields().
     * @param geometry Mesh and finite-element infrastructure.
     * @param domain_parameters Electrode masks, group labels, and global totals.
     * @param bc Boundary conditions for the potential solver.
     * @param cfg Materials, initial values, and cell mode.
     * @param electrode Electrode to initialize (ANODE or CATHODE).
     * @pre The material list is nonempty and matches any particle groups present.
     */
    void Initialize(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        BoundaryConditions& bc, const SimulationConfig& cfg, sim::Electrode electrode);

    /**
     * @brief Evaluate each pair's material-specific chemical potentials.
     *
     * Fields are reset to zero and evaluated only where the pair interface
     * field exceeds 1000.0.
     * @param geometry Provides the number of local vertices.
     * @param avp_pairs Pair interface fields indexed by [j][k], j < k.
     */
    void UpdatePairChemicalPotentials(Initialize_Geometry& geometry,
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs);

    /**
     * @brief Collect aligned field pointers and materials for potential assembly.
     * @param psi Particle-group phase masks in particle order.
     * @param[out] cn_fields Cleared and filled with borrowed concentration pointers.
     * @param[out] psi_fields Cleared and filled with borrowed phase-mask pointers.
     * @param[out] materials Cleared and filled with each group's material.
     * @pre psi contains a valid mask for each particle group.
     */
    void BuildParticleFields(const std::vector<std::unique_ptr<mfem::ParGridFunction>>& psi,
        std::vector<mfem::ParGridFunction*>& cn_fields, std::vector<mfem::ParGridFunction*>& psi_fields, std::vector<sim::MaterialType>& materials) const;

    /**
     * @brief Update each group's exchange current using its material and concentration.
     * @param AvEs Electrode-electrolyte interface fields in particle-group order.
     */
    void UpdateExchangeCurrentDensity(const std::vector<std::unique_ptr<mfem::ParGridFunction>>& AvEs);

    /**
     * @brief Integrate each group's reaction current and sum the results.
     * @param[out] particle_currents Resized and filled with per-group currents.
     * @return Total current over this electrode's particle groups.
     */
    double CalculateElectrodeCurrent(std::vector<double>& particle_currents);

    /**
     * @brief Freeze reaction sources and advance each group's concentration.
     *
     * Copies each current Rxn_gf into Rx_src before the concentration solve.
     * The combined source here is an unweighted sum of these frozen fields.
     * @param weight_pairs Pair coupling weights.
     * @param avp_pairs Pair interface fields.
     * @param ps Phase masks in particle-group order.
     * @param gtPs Global phase totals per group (for example, gtPsA for anodes).
     * @param weightEs Electrode-electrolyte weights in particle-group order.
     * @param[out] total_rxn Reset and filled with the sum of frozen reaction sources.
     */
    void UpdateParticleConcentrations(
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs,
        const std::vector<std::unique_ptr<mfem::ParGridFunction>>& ps, const std::vector<double>& gtPs,
        const std::vector<std::unique_ptr<mfem::ParGridFunction>>& weightEs, mfem::ParGridFunction& total_rxn);

    /**
     * @brief Recompute group reactions and sum their interface-weighted sources.
     * @param[out] total_rxn Reset and filled with the sum of Rxn_gf * WeightEs.
     * @param CnE Electrolyte concentration field.
     * @param phS Shared solid-phase potential field.
     * @param phE Electrolyte potential field.
     * @param AvEs Electrode-electrolyte interface fields in particle-group order.
     * @param WeightEs Electrode-electrolyte weights in particle-group order.
     */
    void UpdateButlerVolmerReactions(mfem::ParGridFunction& total_rxn,
        mfem::ParGridFunction& CnE, mfem::ParGridFunction& phS, mfem::ParGridFunction& phE,
        const std::vector<std::unique_ptr<mfem::ParGridFunction>>& AvEs, const std::vector<std::unique_ptr<mfem::ParGridFunction>>& WeightEs);
};

/**
 * @brief Owns electrode, electrolyte, and combined runtime fields.
 *
 * Full cells initialize both electrodes; half cells initialize only the
 * selected electrode. Solver pointers in the inactive electrode remain null.
 */
struct SimulationState
{
    ElectrodeState anode; ///< Anode state; initialized for full cells and anode half cells.
    ElectrodeState cathode; ///< Cathode state; initialized for full cells and cathode half cells.

    std::unique_ptr<ConcentrationBase> electrolyte_concentration; ///< Electrolyte concentration solver.
    std::unique_ptr<ElectrolytePotential> electrolyte_potential; ///< Electrolyte potential solver.
    std::unique_ptr<mfem::ParGridFunction> CnE_gf; ///< Electrolyte concentration field.
    std::unique_ptr<mfem::ParGridFunction> CnE_gf_psi; ///< Workspace for phase-masked electrolyte concentration.
    std::unique_ptr<mfem::ParGridFunction> phE_gf; ///< Electrolyte potential field.

    std::unique_ptr<Reaction> reaction; ///< Owned reaction helper.
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf; ///< Reaction field; per-group in ParticleState, combined in half-cell SimulationState.
    std::unique_ptr<mfem::ParGridFunction> RxnA_gf; ///< Combined anode reaction source.
    std::unique_ptr<mfem::ParGridFunction> RxnC_gf; ///< Combined cathode reaction source.
    std::unique_ptr<mfem::ParGridFunction> RxnE_gf; ///< Combined electrode reaction source for the full-cell electrolyte.
    std::unique_ptr<mfem::ParGridFunction> CnP_together; ///< Workspace combining particle concentration fields.

    /**
     * @brief Allocate shared fields and initialize the active electrodes.
     *
     * Creates electrolyte solvers and reaction fields, resets both electrode
     * states, and initializes one electrode for HALF mode or both for FULL mode.
     * @param geometry Mesh and finite-element infrastructure.
     * @param domain_parameters Phase masks, particle labels, and global totals.
     * @param bc Prepared boundary conditions.
     * @param cfg Cell mode, active half-cell electrode, materials, and initial values.
     */
    void InitializeFields(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        BoundaryConditions& bc, const SimulationConfig& cfg);
};

#endif // SIMULATION_STATE_HPP
