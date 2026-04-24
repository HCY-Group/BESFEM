#include "../include/SimulationState.hpp"

static double GetInitialValue(
    const std::vector<double>& values,
    int k,
    double fallback)
{
    if (k < static_cast<int>(values.size()))
    {
        return values[k];
    }
    return fallback;
}

static void InitializePairWorkspaces(SimulationState& state, Initialize_Geometry& geometry, int np)
{
    state.mu_pair_a.clear();
    state.mu_pair_b.clear();
    state.sum_pairs.clear();

    state.mu_pair_a.resize(np);
    state.mu_pair_b.resize(np);
    state.sum_pairs.resize(np);

    for (int j = 0; j < np; ++j)
    {
        state.mu_pair_a[j].resize(np);
        state.mu_pair_b[j].resize(np);
        state.sum_pairs[j].resize(np);

        for (int k = 0; k < np; ++k)
        {
            if (j < k)
            {
                state.mu_pair_a[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                state.mu_pair_b[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                state.sum_pairs[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            }
        }
    }

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << "[DEBUG] Initialized pair workspaces for np = "
                  << np << std::endl;
    }
}

static inline double CathodeMuFromC(double c)
{
    return -Constants::Frd * ((1.095 * c * c) - (8.234e-7 * std::exp(14.31 * c)) + (4.692 * std::exp(-0.5389 * c)));
}

void UpdateCathodePairChemicalPotentials(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters)
{
    const int np = static_cast<int>(state.cathode_particles.size());

    for (int j = 0; j < np; ++j)
    {
        for (int k = j + 1; k < np; ++k)
        {
            auto& Cj = *state.cathode_particles[j].Cn_gf;
            auto& Ck = *state.cathode_particles[k].Cn_gf;

            auto& mu_j = *state.mu_pair_a[j][k];
            auto& mu_k = *state.mu_pair_b[j][k];
            auto& pair_if = *domain_parameters.AvP_Pairs[j][k];

            mu_j = 0.0;
            mu_k = 0.0;

            for (int vi = 0; vi < geometry.nV; ++vi)
            {
                if (pair_if(vi) > 1000.0)
                {
                    mu_j(vi) = CathodeMuFromC(Cj(vi));
                    mu_k(vi) = CathodeMuFromC(Ck(vi));
                }
            }
        }
    }
}


static void InitializeAnodeParticles(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    const std::vector<double>& init_values, BoundaryConditions& bc)
{
    const int np = static_cast<int>(domain_parameters.ps.size());
    state.anode_particles.clear();
    state.anode_particles.resize(np);

    if (np == 0)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] No different anode particles defined in the configuration." << std::endl;
        }
        return;
    }

    for (int k = 0; k < np; ++k)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] Creating Anode Particle " << k
                    << " (label = " << domain_parameters.particle_labels[k] << ")"
                    << std::endl;
        }

        auto& p = state.anode_particles[k];
        p.label = domain_parameters.particle_labels[k];

        p.concentration = std::make_unique<CnA>(geometry, domain_parameters);
        p.Cn_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Cn_gf_psi     = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction      = std::make_unique<Reaction>(geometry, domain_parameters);
        p.Rxn_gf        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Rx_src        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.potential     = std::make_unique<PotA>(geometry, domain_parameters, bc);
        p.ph_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        
        p.reaction->Initialize(*p.Rxn_gf, Constants::init_Rxn);

        const double init_cn = GetInitialValue(init_values, k, Constants::init_CnA);

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG]   Initial concentration = " << init_cn << std::endl;
        }

        p.concentration->SetupField(*p.Cn_gf, init_cn, *domain_parameters.ps[k], domain_parameters.gtPs[k]);
        p.potential->SetupField(*p.ph_gf, Constants::init_BvA, *domain_parameters.ps[k]);
    }
}

static void InitializeCathodeParticles(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    const std::vector<double>& init_values, BoundaryConditions& bc)
{
    const int np = static_cast<int>(domain_parameters.ps.size());
    state.cathode_particles.clear();
    state.cathode_particles.resize(np);

    if (np == 0)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] No different cathode particles defined in the configuration." << std::endl;
        }
        return;
    }

    for (int k = 0; k < np; ++k)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] Creating Cathode Particle " << k
                    << " (label = " << domain_parameters.particle_labels[k] << ")"
                    << std::endl;
        }

        auto& p = state.cathode_particles[k];
        p.label = domain_parameters.particle_labels[k];

        p.concentration = std::make_unique<CnC>(geometry, domain_parameters);
        p.Cn_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Cn_gf_psi     = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction      = std::make_unique<Reaction>(geometry, domain_parameters);
        p.Rxn_gf        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Rx_src        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.potential     = std::make_unique<PotC>(geometry, domain_parameters, bc);
        p.ph_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction->Initialize(*p.Rxn_gf, Constants::init_Rxn);

        const double init_cn = GetInitialValue(init_values, k, Constants::init_CnC);

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG]   Initial concentration = " << init_cn << std::endl;
        }

        p.concentration->SetupField(*p.Cn_gf, init_cn, *domain_parameters.ps[k], domain_parameters.gtPs[k]);
        p.potential->SetupField(*p.ph_gf, Constants::init_BvC, *domain_parameters.ps[k]);
    }
}

void InitializeFields(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    BoundaryConditions& bc, const SimulationConfig& cfg)
{
    state.CnP_together = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.CnE_gf_psi   = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

    state.electrolyte_concentration = std::make_unique<CnE>(geometry, domain_parameters, bc, cfg.mode);
    state.CnE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.electrolyte_concentration->SetupField(*state.CnE_gf, Constants::init_CnE, *domain_parameters.pse, domain_parameters.gtPse);

    state.electrolyte_potential = std::make_unique<PotE>(geometry, domain_parameters, bc, cfg.mode);
    state.phE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.electrolyte_potential->SetupField(*state.phE_gf, Constants::init_BvE, *domain_parameters.pse);

    state.reaction = std::make_unique<Reaction>(geometry, domain_parameters);
    state.Rxn_gf   = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.reaction->Initialize(*state.Rxn_gf, Constants::init_Rxn);

    if (cfg.mode == sim::CellMode::HALF)
    {
        if (cfg.half_electrode == sim::Electrode::ANODE)
        {
            const int np = static_cast<int>(domain_parameters.ps.size());

            state.CnA_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

            state.anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
            state.CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.anode_concentration->SetupField(*state.CnA_gf, Constants::init_CnA, *domain_parameters.psi, domain_parameters.gtPsi);

            state.anode_potential = std::make_unique<PotA>(geometry, domain_parameters, bc);
            state.phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.anode_potential->SetupField(*state.phA_gf, Constants::init_BvA, *domain_parameters.psi);

            InitializeAnodeParticles(state, geometry, domain_parameters, cfg.init_anode_particles, bc);
            InitializePairWorkspaces(state, geometry, static_cast<int>(state.anode_particles.size()));

            state.anode_out.clear();
            state.anode_out.resize(np);

            for (int k = 0; k < np; ++k)
            {
                state.anode_out[k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            }

        }
        else
        {
            const int np = static_cast<int>(domain_parameters.ps.size());

            state.CnC_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

            state.cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
            state.CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.cathode_concentration->SetupField(*state.CnC_gf, Constants::init_CnC, *domain_parameters.psi, domain_parameters.gtPsi);

            state.cathode_potential = std::make_unique<PotC>(geometry, domain_parameters, bc);
            state.phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.cathode_potential->SetupField(*state.phC_gf, Constants::init_BvC, *domain_parameters.psi);

            InitializeCathodeParticles(state, geometry, domain_parameters, cfg.init_cathode_particles, bc);
            InitializePairWorkspaces(state, geometry,static_cast<int>(state.cathode_particles.size()));

            state.cathode_out.clear();
            state.cathode_out.resize(np);

            for (int k = 0; k < np; ++k)
            {
                state.cathode_out[k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            }
        }
    }
    else
    {
        state.CnA_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.CnC_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        state.anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
        state.CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.anode_concentration->SetupField(*state.CnA_gf, Constants::init_CnA, *domain_parameters.psA, domain_parameters.gtPsA);

        state.anode_potential = std::make_unique<PotA>(geometry, domain_parameters, bc);
        state.phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.anode_potential->SetupField(*state.phA_gf, Constants::init_BvA, *domain_parameters.psA);

        state.cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
        state.CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.cathode_concentration->SetupField(*state.CnC_gf, Constants::init_CnC, *domain_parameters.psC, domain_parameters.gtPsC);

        state.cathode_potential = std::make_unique<PotC>(geometry, domain_parameters, bc);
        state.phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.cathode_potential->SetupField(*state.phC_gf, Constants::init_BvC, *domain_parameters.psC);
    }

    if (mfem::Mpi::WorldRank() == 0 && (state.anode_particles.size() > 0 || state.cathode_particles.size() > 0))
    {
        std::cout << "[DEBUG] Finished InitializeFields()" << std::endl;

        std::cout << "    Anode particles:   "
                << state.anode_particles.size() << std::endl;

        std::cout << "    Cathode particles: "
                << state.cathode_particles.size() << std::endl;
    }
}

void Pairs(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters, int j, std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np)
{                        
    for (int k = 0; k < np; ++k)
    {
        if (j == k) { continue; }

        const int a = std::min(j, k);
        const int b = std::max(j, k);

        ConcentrationBase::PairCoupling pair;
        pair.sum_part = state.sum_pairs[a][b].get();
        pair.weight   = domain_parameters.WeightPairs[a][b].get();
        pair.grad_psi = domain_parameters.AvP_Pairs[a][b].get();

        if (j < k)
        {
            pair.mu_self = state.mu_pair_a[a][b].get();
            pair.mu_nbr  = state.mu_pair_b[a][b].get();
        }
        else
        {
            pair.mu_self = state.mu_pair_b[a][b].get();
            pair.mu_nbr  = state.mu_pair_a[a][b].get();
        }

        pair_terms.push_back(pair);

        if (mfem::Mpi::WorldRank() == 0 && t == 1)
        {
            std::cout << "[DEBUG] Pair (j,k) = (" << j << "," << k << ")"
                    << " | stored as (a,b) = (" << a << "," << b << ")"; 
            std::cout << std::endl;
        }

    }
}
