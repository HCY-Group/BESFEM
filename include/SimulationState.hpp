/** Runtime fields and the methods that initialize and update them. */
#ifndef SIMULATION_STATE_HPP
#define SIMULATION_STATE_HPP

#include "mfem.hpp"
#include "BESFEM_All.hpp"

#include <memory>
#include <vector>

// One particle group's concentration and reaction solvers and fields.
struct ParticleState
{
    int label = -1;
    sim::MaterialType material = sim::MaterialType::Graphite;
    std::unique_ptr<ConcentrationBase> concentration;
    std::unique_ptr<mfem::ParGridFunction> Cn_gf;
    std::unique_ptr<mfem::ParGridFunction> Cn_gf_psi;
    std::unique_ptr<Reaction> reaction;
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf;
    std::unique_ptr<mfem::ParGridFunction> Rx_src;

    void Initialize(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        const SimulationConfig& cfg, sim::MaterialType particle_material, int particle_label,
        double init_cn, mfem::ParGridFunction& particle_field, double particle_total);
};

// Fields shared by each pair of particle groups within one electrode.
struct PairWorkspaces
{
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_a;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_b;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> sum_pairs;

    void Initialize(Initialize_Geometry& geometry, int np, const char* electrode_name);

    void BuildPairTerms(const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs, int j,
        std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np) const;
};

// Anode and cathode use the same methods, but own separate fields and solvers.
struct ElectrodeState
{
    std::vector<ParticleState> particles;
    PairWorkspaces pairs;
    std::unique_ptr<ElectrodePotential> potential;
    std::unique_ptr<mfem::ParGridFunction> ph_gf;

    void Initialize(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        BoundaryConditions& bc, const SimulationConfig& cfg, sim::Electrode electrode);

    void UpdatePairChemicalPotentials(Initialize_Geometry& geometry,
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs);

    void BuildParticleFields(const std::vector<std::unique_ptr<mfem::ParGridFunction>>& psi,
        std::vector<mfem::ParGridFunction*>& cn_fields, std::vector<mfem::ParGridFunction*>& psi_fields, std::vector<sim::MaterialType>& materials) const;

    void UpdateExchangeCurrentDensity(const std::vector<std::unique_ptr<mfem::ParGridFunction>>& AvEs);

    double CalculateElectrodeCurrent(std::vector<double>& particle_currents);

    void UpdateParticleConcentrations(
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
        const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs,
        const std::vector<std::unique_ptr<mfem::ParGridFunction>>& ps, const std::vector<double>& gtPs,
        const std::vector<std::unique_ptr<mfem::ParGridFunction>>& weightEs, mfem::ParGridFunction& total_rxn);

    void UpdateButlerVolmerReactions(mfem::ParGridFunction& total_rxn,
        mfem::ParGridFunction& CnE, mfem::ParGridFunction& phS, mfem::ParGridFunction& phE,
        const std::vector<std::unique_ptr<mfem::ParGridFunction>>& AvEs, const std::vector<std::unique_ptr<mfem::ParGridFunction>>& WeightEs);
};

// Owns both electrodes, the electrolyte, and the combined reaction fields.
struct SimulationState
{
    ElectrodeState anode;
    ElectrodeState cathode;

    std::unique_ptr<ConcentrationBase> electrolyte_concentration;
    std::unique_ptr<ElectrolytePotential> electrolyte_potential;
    std::unique_ptr<mfem::ParGridFunction> CnE_gf;
    std::unique_ptr<mfem::ParGridFunction> CnE_gf_psi;
    std::unique_ptr<mfem::ParGridFunction> phE_gf;

    std::unique_ptr<Reaction> reaction;
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf;
    std::unique_ptr<mfem::ParGridFunction> RxnA_gf;
    std::unique_ptr<mfem::ParGridFunction> RxnC_gf;
    std::unique_ptr<mfem::ParGridFunction> RxnE_gf;
    std::unique_ptr<mfem::ParGridFunction> CnP_together;

    void InitializeFields(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        BoundaryConditions& bc, const SimulationConfig& cfg);
};

#endif // SIMULATION_STATE_HPP
