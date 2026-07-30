/**
 * @file SimulationState.hpp
 * @brief Defines runtime state containers for BESFEM simulations.
 */

#ifndef SIMULATION_STATE_HPP
#define SIMULATION_STATE_HPP

#include "mfem.hpp"
#include "BESFEM_All.hpp"

#include <memory>
#include <vector>

/**
 * @struct AnodeParticleState
 * @brief Stores solver objects and fields for one anode particle group.
 *
 * AnodeParticleState collects the concentration solver, reaction solver,
 * potential solver, and associated MFEM grid functions for a single labeled
 * anode particle or particle group.
 */
struct AnodeParticleState
{
    int label = -1; ///< Integer label identifying this particle in the input geometry.

    sim::MaterialType material = sim::MaterialType::Graphite; ///< Material assigned to this anode particle.

    std::unique_ptr<ConcentrationBase> concentration; ///< Concentration solver for this particle.
    std::unique_ptr<mfem::ParGridFunction> Cn_gf; ///< Particle concentration field.
    std::unique_ptr<mfem::ParGridFunction> Cn_gf_psi; ///< Concentration field masked by the particle phase field.

    std::unique_ptr<Reaction> reaction; ///< Reaction model for this particle.
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf; ///< Butler--Volmer reaction field.
    std::unique_ptr<mfem::ParGridFunction> Rx_src; ///< Reaction source term used by the concentration update.

    std::unique_ptr<PotentialBase> potential; ///< Solid-phase potential solver for this particle.
    std::unique_ptr<mfem::ParGridFunction> ph_gf; ///< Solid-phase potential field.
};

/**
 * @struct CathodeParticleState
 * @brief Stores solver objects and fields for one cathode particle group.
 *
 * CathodeParticleState collects the concentration solver, reaction solver,
 * potential solver, and associated MFEM grid functions for a single labeled
 * cathode particle or particle group.
 */
struct CathodeParticleState
{
    int label = -1; ///< Integer label identifying this particle in the input geometry.

    sim::MaterialType material = sim::MaterialType::NMC; ///< Material assigned to this cathode particle.

    std::unique_ptr<ConcentrationBase> concentration; ///< Concentration solver for this particle.
    std::unique_ptr<mfem::ParGridFunction> Cn_gf; ///< Particle concentration field.
    std::unique_ptr<mfem::ParGridFunction> Cn_gf_psi; ///< Concentration field masked by the particle phase field.

    std::unique_ptr<Reaction> reaction; ///< Reaction model for this particle.
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf; ///< Butler--Volmer reaction field.
    std::unique_ptr<mfem::ParGridFunction> Rx_src; ///< Reaction source term used by the concentration update.

    std::unique_ptr<PotentialBase> potential; ///< Solid-phase potential solver for this particle.
    std::unique_ptr<mfem::ParGridFunction> ph_gf; ///< Solid-phase potential field.
};

struct PairWorkspaces
{
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_a;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_b;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> sum_pairs;
};

/**
 * @struct SimulationState
 * @brief Owns the runtime solvers and grid functions for a BESFEM simulation.
 *
 * SimulationState is the main container for all dynamically allocated solver
 * objects and MFEM fields used during a run. It stores single-electrode fields,
 * electrolyte fields, particle-resolved fields, reaction fields, pairwise
 * chemical-potential workspaces, and output workspaces.
 */
struct SimulationState
{
    std::unique_ptr<ConcentrationBase> anode_concentration; ///< Anode concentration solver.
    // std::unique_ptr<PotentialBase> anode_potential; ///< Anode solid-potential solver.
    std::shared_ptr<ElectrodePotential> anode_potential; ///< Anode solid-potential solver.
    std::unique_ptr<mfem::ParGridFunction> CnA_gf; ///< Anode concentration field.
    std::unique_ptr<mfem::ParGridFunction> CnA_gf_psi; ///< Masked anode concentration field.
    std::unique_ptr<mfem::ParGridFunction> phA_gf; ///< Anode solid-potential field.

    std::unique_ptr<ConcentrationBase> cathode_concentration; ///< Cathode concentration solver.
    // std::unique_ptr<PotentialBase> cathode_potential; ///< Cathode solid-potential solver.
    std::shared_ptr<ElectrodePotential> cathode_potential; ///< Cathode solid-potential solver.
    std::unique_ptr<mfem::ParGridFunction> CnC_gf; ///< Cathode concentration field.
    std::unique_ptr<mfem::ParGridFunction> CnC_gf_psi; ///< Masked cathode concentration field.
    std::unique_ptr<mfem::ParGridFunction> phC_gf; ///< Cathode solid-potential field.

    std::unique_ptr<ConcentrationBase> electrolyte_concentration; ///< Electrolyte concentration solver.
    // std::unique_ptr<ElectrodePotential> electrolyte_potential; ///< Electrolyte potential solver.
    std::unique_ptr<ElectrolytePotential> electrolyte_potential; ///< Electrolyte potential solver.
    std::unique_ptr<mfem::ParGridFunction> CnE_gf; ///< Electrolyte concentration field.
    std::unique_ptr<mfem::ParGridFunction> CnE_gf_psi; ///< Masked electrolyte concentration field.
    std::unique_ptr<mfem::ParGridFunction> phE_gf; ///< Electrolyte potential field.

    std::unique_ptr<Reaction> reaction; ///< Global reaction model for non-particle-resolved simulations.
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf; ///< Global reaction field.
    std::unique_ptr<mfem::ParGridFunction> RxnA_gf; ///< Anode reaction field.
    std::unique_ptr<mfem::ParGridFunction> RxnC_gf; ///< Cathode reaction field.
    std::unique_ptr<mfem::ParGridFunction> RxnE_gf; ///< Electrolyte reaction field.

    std::unique_ptr<mfem::ParGridFunction> CnP_together; ///< Combined particle concentration field.

    std::vector<AnodeParticleState> anode_particles; ///< Per-particle anode state objects.
    std::vector<CathodeParticleState> cathode_particles; ///< Per-particle cathode state objects.

    PairWorkspaces anode_pairs;
    PairWorkspaces cathode_pairs;

    // std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_a; ///< Pairwise chemical-potential fields for particle A.
    // std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_b; ///< Pairwise chemical-potential fields for particle B.
    // std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> sum_pairs; ///< Pairwise accumulated coupling fields.

    std::vector<std::unique_ptr<mfem::ParGridFunction>> cathode_out; ///< Output workspaces for cathode concentration fields.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> anode_out; ///< Output workspaces for anode concentration fields.
};


/**
 * @brief Allocate and initialize all simulation fields and solver objects.
 *
 * Creates concentration, potential, and reaction solvers according to the
 * simulation configuration, allocates their associated grid functions, and
 * initializes field values before the timestep loop begins.
 *
 * @param state Simulation state object to populate.
 * @param geometry Geometry and finite-element-space handler.
 * @param domain_parameters Domain fields and global quantities.
 * @param bc Boundary-condition handler.
 * @param cfg Simulation configuration.
 */
void InitializeFields(SimulationState& state,
                      Initialize_Geometry& geometry,
                      Domain_Parameters& domain_parameters,
                      BoundaryConditions& bc,
                      const SimulationConfig& cfg);

/**
 * @brief Update pairwise cathode chemical-potential fields.
 *
 * Computes chemical-potential workspaces used for cathode particle-particle
 * coupling terms.
 *
 * @param state Simulation state containing cathode particle fields.
 * @param geometry Geometry and finite-element-space handler.
 * @param avp_pairs Pairwise chemical potential fields.
 */
void UpdateCathodePairChemicalPotentials(SimulationState& state,
                                         Initialize_Geometry& geometry,
                                         const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs);

/**
 * @brief Update pairwise anode chemical-potential fields.
 *
 * Computes chemical-potential workspaces used for anode particle-particle
 * coupling terms.
 *
 * @param state Simulation state containing anode particle fields.
 * @param geometry Geometry and finite-element-space handler.
 * @param avp_pairs Pairwise chemical potential fields.
 */
void UpdateAnodePairChemicalPotentials(SimulationState& state,
                                       Initialize_Geometry& geometry,
                                       const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs);

// /**
//  * @brief Build pairwise coupling terms for one particle.
//  *
//  * Populates the pair_terms vector with the interface weights, chemical
//  * potentials, and accumulated pair fields needed by the concentration solver
//  * for particle-particle coupling.
//  *
//  * @param state Simulation state containing pairwise workspaces.
//  * @param geometry Geometry and finite-element-space handler.
//  * @param domain_parameters Domain fields and pairwise coupling masks.
//  * @param j Particle index for which pair terms are assembled.
//  * @param pair_terms Output vector of pairwise coupling structures.
//  * @param np Number of particles in the current electrode group.
//  * @param t Current timestep index.
//  */
// void Pairs(SimulationState& state,
//            Initialize_Geometry& geometry,
//            Domain_Parameters& domain_parameters,
//            int j,
//            std::vector<ConcentrationBase::PairCoupling>& pair_terms,
//            int np,
//            int t);

void Pairs(PairWorkspaces& workspace, const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
    const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs,
    int j, std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np, int t);

#endif // SIMULATION_STATE_HPP
