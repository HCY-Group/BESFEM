#ifndef SIMULATION_STATE_HPP
#define SIMULATION_STATE_HPP

#include "mfem.hpp"
#include "BESFEM_All.hpp"
#include <memory>
#include <vector>

struct AnodeParticleState
{
    int label = -1;

    std::unique_ptr<CnA> concentration;
    std::unique_ptr<mfem::ParGridFunction> Cn_gf;
    std::unique_ptr<mfem::ParGridFunction> Cn_gf_psi;

    std::unique_ptr<Reaction> reaction;
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf;
    std::unique_ptr<mfem::ParGridFunction> Rx_src;

    std::unique_ptr<PotA> potential;
    std::unique_ptr<mfem::ParGridFunction> ph_gf;
};

struct CathodeParticleState
{
    int label = -1;

    std::unique_ptr<CnC> concentration;
    std::unique_ptr<mfem::ParGridFunction> Cn_gf;
    std::unique_ptr<mfem::ParGridFunction> Cn_gf_psi;

    std::unique_ptr<Reaction> reaction;
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf;
    std::unique_ptr<mfem::ParGridFunction> Rx_src;

    std::unique_ptr<PotC> potential;
    std::unique_ptr<mfem::ParGridFunction> ph_gf;
};

struct SimulationState
{
    std::unique_ptr<CnA> anode_concentration;
    std::unique_ptr<PotA> anode_potential;
    std::unique_ptr<mfem::ParGridFunction> CnA_gf, CnA_gf_psi, phA_gf;

    std::unique_ptr<CnC> cathode_concentration;
    std::unique_ptr<PotC> cathode_potential;
    std::unique_ptr<mfem::ParGridFunction> CnC_gf, CnC_gf_psi, phC_gf;

    std::unique_ptr<CnE> electrolyte_concentration;
    std::unique_ptr<PotE> electrolyte_potential;
    std::unique_ptr<mfem::ParGridFunction> CnE_gf, CnE_gf_psi, phE_gf;

    std::unique_ptr<Reaction> reaction;
    std::unique_ptr<mfem::ParGridFunction> Rxn_gf;

    std::unique_ptr<mfem::ParGridFunction> CnP_together;

    std::vector<AnodeParticleState> anode_particles;
    std::vector<CathodeParticleState> cathode_particles;

    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_a;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> mu_pair_b;
    std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> sum_pairs;

    std::vector<std::unique_ptr<mfem::ParGridFunction>> cathode_out;
    std::vector<std::unique_ptr<mfem::ParGridFunction>> anode_out;
};

void InitializeFields(
    SimulationState& state,
    Initialize_Geometry& geometry,
    Domain_Parameters& domain_parameters,
    BoundaryConditions& bc,
    const SimulationConfig& cfg);

void UpdateCathodePairChemicalPotentials(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters);

void Pairs(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters, int j, std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np);

#endif