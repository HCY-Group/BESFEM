#include "../include/SimulationState.hpp"
#include "../include/MaterialProperties.hpp"

void PairWorkspaces::Initialize(Initialize_Geometry& geometry, int np, const char* electrode_name)
{
    mu_pair_a.clear();
    mu_pair_b.clear();
    sum_pairs.clear();

    mu_pair_a.resize(np);
    mu_pair_b.resize(np);
    sum_pairs.resize(np);

    for (int j = 0; j < np; ++j)
    {
        mu_pair_a[j].resize(np);
        mu_pair_b[j].resize(np);
        sum_pairs[j].resize(np);

        for (int k = j + 1; k < np; ++k)
        {
            mu_pair_a[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            mu_pair_b[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            sum_pairs[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

            *mu_pair_a[j][k] = 0.0;
            *mu_pair_b[j][k] = 0.0;
            *sum_pairs[j][k] = 0.0;
        }
    }

    if (mfem::Mpi::WorldRank() == 0)
    {
        const int number_of_pairs = np * (np - 1) / 2;

        std::cout << "[DEBUG] Initialized " << electrode_name << " pair workspaces for np = " << np << " (" << number_of_pairs << " pairs)" << std::endl;
    }
}

void ParticleState::Initialize(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    const SimulationConfig& cfg, sim::MaterialType particle_material, int particle_label,
    double init_cn, mfem::ParGridFunction& particle_field, double particle_total)
{
    label = particle_label;
    material = particle_material;

    switch (material)
    {
        case sim::MaterialType::Graphite:
        case sim::MaterialType::LFP:
            concentration = std::make_unique<ElectrodeCahnHilliard>(geometry, domain_parameters, material, cfg);
            break;
        case sim::MaterialType::Carbon:
        case sim::MaterialType::Silicon:
        case sim::MaterialType::NMC:
            concentration = std::make_unique<ElectrodeDiffusion>(geometry, domain_parameters, material, cfg);
            break;
        default:
            mfem::mfem_error("Unsupported electrode material physics.");
    }

    Cn_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    Cn_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    reaction = std::make_unique<Reaction>(geometry, domain_parameters, cfg);
    Rxn_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    Rx_src = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

    reaction->Initialize(*Rxn_gf, Constants::init_Rxn);
    concentration->SetupField(*Cn_gf, init_cn, particle_field, particle_total);
}

void ElectrodeState::Initialize(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    BoundaryConditions& bc, const SimulationConfig& cfg, sim::Electrode electrode)
{
    const bool is_anode = electrode == sim::Electrode::ANODE;
    const bool is_half = cfg.mode == sim::CellMode::HALF;
    const char* electrode_name = is_anode ? "anode" : "cathode";

    const auto& materials = is_anode ? cfg.anode_materials : cfg.cathode_materials;
    const auto& init_values = is_anode ? cfg.init_anode_particles : cfg.init_cathode_particles;
    const double init_cn = is_anode ? cfg.init_CnA : cfg.init_CnC;
    const double init_bv = is_anode ? cfg.init_BvA : cfg.init_BvC;

    const auto& particle_fields = is_half ? domain_parameters.ps : (is_anode ? domain_parameters.psA : domain_parameters.psC);
    const auto& particle_totals = is_half ? domain_parameters.gtPs : (is_anode ? domain_parameters.gtPsA : domain_parameters.gtPsC);
    const auto& particle_labels = is_half ? domain_parameters.particle_labels : (is_anode ? domain_parameters.anode_particle_labels : domain_parameters.cathode_particle_labels);
    auto& psi = is_half ? domain_parameters.psi : (is_anode ? domain_parameters.psiA : domain_parameters.psiC);
    const int np = static_cast<int>(particle_fields.size());

    MFEM_VERIFY(!materials.empty(), "An electrode material must be specified before initializing potential.");
    MFEM_VERIFY(np == 0 || materials.size() == particle_fields.size(), "Provide one material for each particle group.");
    if (!is_anode && np > 0)
    {
        MFEM_VERIFY(init_values.size() == particle_fields.size(), "Provide one initial concentration for each cathode particle group.");
    }

    potential = std::make_unique<ElectrodePotential>(geometry, domain_parameters, bc, electrode, materials.front(), cfg);
    ph_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    potential->SetupField(*ph_gf, init_bv, *psi);

    particles.clear();
    particles.resize(np);
    for (int k = 0; k < np; ++k)
    {
        const sim::MaterialType material = materials[k];
        if (is_anode)
        {
            MFEM_VERIFY(material == sim::MaterialType::Graphite || material == sim::MaterialType::Carbon || material == sim::MaterialType::Silicon, "Unsupported anode material.");
        }
        else
        {
            MFEM_VERIFY(material == sim::MaterialType::NMC || material == sim::MaterialType::LFP, "Unsupported cathode material.");
        }

        double particle_init_cn = init_cn;
        if (k < static_cast<int>(init_values.size()))
        {
            particle_init_cn = init_values[k];
        }
        particles[k].Initialize(geometry, domain_parameters, cfg, material, particle_labels[k], particle_init_cn, *particle_fields[k], particle_totals[k]);
    }
    pairs.Initialize(geometry, np, electrode_name);
}

void SimulationState::InitializeFields(Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    BoundaryConditions& bc, const SimulationConfig& cfg)
{
    CnP_together = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    CnE_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

    electrolyte_concentration = std::make_unique<ElectrolyteDiffusion>(geometry, domain_parameters, bc, cfg.mode, sim::MaterialType::Electrolyte, cfg);
    CnE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    electrolyte_concentration->SetupField(*CnE_gf, cfg.init_CnE, *domain_parameters.pse, domain_parameters.gtPse);

    electrolyte_potential = std::make_unique<ElectrolytePotential>(geometry, domain_parameters, bc, cfg.mode, cfg);
    phE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    electrolyte_potential->SetupField(*phE_gf, cfg.init_BvE, *domain_parameters.pse);

    reaction = std::make_unique<Reaction>(geometry, domain_parameters, cfg);
    Rxn_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    RxnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    RxnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    RxnE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

    reaction->Initialize(*Rxn_gf, Constants::init_Rxn);
    reaction->Initialize(*RxnA_gf, Constants::init_Rxn);
    reaction->Initialize(*RxnC_gf, Constants::init_Rxn);
    reaction->Initialize(*RxnE_gf, Constants::init_Rxn);

    // A half cell initializes one electrode; a full cell initializes both.
    anode = ElectrodeState{};
    cathode = ElectrodeState{};
    if (cfg.mode == sim::CellMode::FULL || cfg.half_electrode == sim::Electrode::ANODE)
    {
        anode.Initialize(geometry, domain_parameters, bc, cfg, sim::Electrode::ANODE);
    }
    if (cfg.mode == sim::CellMode::FULL || cfg.half_electrode == sim::Electrode::CATHODE)
    {
        cathode.Initialize(geometry, domain_parameters, bc, cfg, sim::Electrode::CATHODE);
    }
}

void ElectrodeState::UpdatePairChemicalPotentials(Initialize_Geometry& geometry,
    const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs)
{
    const int np = static_cast<int>(particles.size());

    for (int j = 0; j < np; ++j)
    {
        for (int k = j + 1; k < np; ++k)
        {
            auto& Cj = *particles[j].Cn_gf;
            auto& Ck = *particles[k].Cn_gf;

            const auto mat_j = particles[j].material;
            const auto mat_k = particles[k].material;

            auto& mu_j = *pairs.mu_pair_a[j][k];
            auto& mu_k = *pairs.mu_pair_b[j][k];
            auto& AvP_pair = *avp_pairs[j][k];

            mu_j = 0.0;
            mu_k = 0.0;

            for (int vi = 0; vi < geometry.nV; ++vi)
            {
                if (AvP_pair(vi) > 1000.0)
                {
                    mu_j(vi) = MaterialProperties::ChemicalPotential(mat_j, Cj(vi));
                    mu_k(vi) = MaterialProperties::ChemicalPotential(mat_k, Ck(vi));
                }
            }
        }
    }
}

void PairWorkspaces::BuildPairTerms(const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
    const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs, int j,
    std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np) const
{
    pair_terms.clear();

    for (int k = 0; k < np; ++k)
    {
        if (j == k)
        {
            continue;
        }

        const int a = std::min(j, k);
        const int b = std::max(j, k);

        MFEM_VERIFY(sum_pairs[a][b], "Pair sum workspace is null.");
        MFEM_VERIFY(mu_pair_a[a][b], "Pair chemical-potential workspace A is null.");
        MFEM_VERIFY(mu_pair_b[a][b], "Pair chemical-potential workspace B is null.");
        MFEM_VERIFY(weight_pairs[a][b], "Pair weight field is null.");
        MFEM_VERIFY(avp_pairs[a][b], "Pair interface field is null.");

        ConcentrationBase::PairCoupling pair;

        pair.sum_part = sum_pairs[a][b].get();
        pair.weight = weight_pairs[a][b].get();
        pair.grad_psi = avp_pairs[a][b].get();

        if (j < k)
        {
            pair.mu_self = mu_pair_a[a][b].get();
            pair.mu_nbr = mu_pair_b[a][b].get();
        }
        else
        {
            pair.mu_self = mu_pair_b[a][b].get();
            pair.mu_nbr = mu_pair_a[a][b].get();
        }

        pair_terms.push_back(pair);
    }
}

void ElectrodeState::BuildParticleFields(const std::vector<std::unique_ptr<mfem::ParGridFunction>>& psi,
    std::vector<mfem::ParGridFunction*>& cn_fields, std::vector<mfem::ParGridFunction*>& psi_fields, std::vector<sim::MaterialType>& materials) const
{
    const int np = static_cast<int>(particles.size());

    cn_fields.clear();
    psi_fields.clear();
    materials.clear();

    cn_fields.reserve(np);
    psi_fields.reserve(np);
    materials.reserve(np);

    for (int j = 0; j < np; ++j)
    {
        cn_fields.push_back(particles[j].Cn_gf.get());
        psi_fields.push_back(psi[j].get());
        materials.push_back(particles[j].material);
    }
}

void ElectrodeState::UpdateExchangeCurrentDensity(const std::vector<std::unique_ptr<mfem::ParGridFunction>>& AvEs)
{
    const int np = static_cast<int>(particles.size());

    for (int j = 0; j < np; ++j)
    {
        particles[j].reaction->ExchangeCurrentDensity(*particles[j].Cn_gf, *AvEs[j], particles[j].material);
    }
}

double ElectrodeState::CalculateElectrodeCurrent(std::vector<double>& particle_currents)
{
    const int np = static_cast<int>(particles.size());

    particle_currents.assign(np, 0.0);

    double total_current = 0.0;

    for (int j = 0; j < np; ++j)
    {
        particles[j].reaction->TotalReactionCurrent(*particles[j].Rxn_gf, particle_currents[j]);
        total_current += particle_currents[j];
    }

    return total_current;
}

void ElectrodeState::UpdateParticleConcentrations(
    const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
    const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs,
    const std::vector<std::unique_ptr<mfem::ParGridFunction>>& ps, const std::vector<double>& gtPs,
    const std::vector<std::unique_ptr<mfem::ParGridFunction>>& weightEs, mfem::ParGridFunction& total_rxn)
{
    const int np = static_cast<int>(particles.size());

    total_rxn = 0.0;

    for (int j = 0; j < np; ++j)
    {
        // Freeze the reaction from the previous converged timestep.
        *particles[j].Rx_src = *particles[j].Rxn_gf;

        total_rxn += *particles[j].Rx_src;
        std::vector<ConcentrationBase::PairCoupling> pair_terms;

        pairs.BuildPairTerms(weight_pairs, avp_pairs, j, pair_terms, np);
        particles[j].concentration->UpdateConcentration(*particles[j].Rx_src, *particles[j].Cn_gf, *ps[j], gtPs[j], *weightEs[j], pair_terms);
    }
}

void ElectrodeState::UpdateButlerVolmerReactions(mfem::ParGridFunction& total_rxn,
    mfem::ParGridFunction& CnE, mfem::ParGridFunction& phS, mfem::ParGridFunction& phE,
    const std::vector<std::unique_ptr<mfem::ParGridFunction>>& AvEs, const std::vector<std::unique_ptr<mfem::ParGridFunction>>& WeightEs)
{
    total_rxn = 0.0;

    const int np = static_cast<int>(particles.size());

    for (int j = 0; j < np; ++j)
    {
        particles[j].reaction->ButlerVolmer(*particles[j].Rxn_gf, *particles[j].Cn_gf, CnE, phS, phE, *AvEs[j]);

        mfem::ParGridFunction weighted_rxn(particles[j].Rxn_gf->ParFESpace());

        weighted_rxn = *particles[j].Rxn_gf;
        weighted_rxn *= *WeightEs[j];

        total_rxn += weighted_rxn;
    }
}
