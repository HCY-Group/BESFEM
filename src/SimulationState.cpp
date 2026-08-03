#include "../include/SimulationState.hpp"
#include "../include/MaterialProperties.hpp"

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

static void InitializePairWorkspaces(
    PairWorkspaces& workspace,
    Initialize_Geometry& geometry,
    int np,
    const char* electrode_name)
{
    workspace.mu_pair_a.clear();
    workspace.mu_pair_b.clear();
    workspace.sum_pairs.clear();

    workspace.mu_pair_a.resize(np);
    workspace.mu_pair_b.resize(np);
    workspace.sum_pairs.resize(np);

    for (int j = 0; j < np; ++j)
    {
        workspace.mu_pair_a[j].resize(np);
        workspace.mu_pair_b[j].resize(np);
        workspace.sum_pairs[j].resize(np);

        for (int k = j + 1; k < np; ++k)
        {
            workspace.mu_pair_a[j][k] =
                std::make_unique<mfem::ParGridFunction>(
                    geometry.parfespace.get());

            workspace.mu_pair_b[j][k] =
                std::make_unique<mfem::ParGridFunction>(
                    geometry.parfespace.get());

            workspace.sum_pairs[j][k] =
                std::make_unique<mfem::ParGridFunction>(
                    geometry.parfespace.get());

            *workspace.mu_pair_a[j][k] = 0.0;
            *workspace.mu_pair_b[j][k] = 0.0;
            *workspace.sum_pairs[j][k] = 0.0;
        }
    }

    if (mfem::Mpi::WorldRank() == 0)
    {
        const int number_of_pairs = np * (np - 1) / 2;

        std::cout
            << "[DEBUG] Initialized "
            << electrode_name
            << " pair workspaces for np = "
            << np
            << " (" << number_of_pairs << " pairs)"
            << std::endl;
    }
}

void UpdateCathodePairChemicalPotentials(SimulationState& state, Initialize_Geometry& geometry, const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs)
{
    const int np = static_cast<int>(state.cathode_particles.size());
    // InitializePairWorkspaces(state.cathode_pairs, geometry, static_cast<int>(state.cathode_particles.size()), "cathode");

    for (int j = 0; j < np; ++j)
    {
        for (int k = j + 1; k < np; ++k)
        {
            auto& Cj = *state.cathode_particles[j].Cn_gf;
            auto& Ck = *state.cathode_particles[k].Cn_gf;

            const auto mat_j = state.cathode_particles[j].material;
            const auto mat_k = state.cathode_particles[k].material;

            auto& mu_j = *state.cathode_pairs.mu_pair_a[j][k];
            auto& mu_k = *state.cathode_pairs.mu_pair_b[j][k];
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

void UpdateAnodePairChemicalPotentials(SimulationState& state, Initialize_Geometry& geometry, const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs)
{
    const int np = static_cast<int>(state.anode_particles.size());
    // InitializePairWorkspaces(state.anode_pairs, geometry, static_cast<int>(state.anode_particles.size()), "anode");

    for (int j = 0; j < np; ++j)
    {
        for (int k = j + 1; k < np; ++k)
        {
            auto& Cj = *state.anode_particles[j].Cn_gf;
            auto& Ck = *state.anode_particles[k].Cn_gf;

            const auto mat_j = state.anode_particles[j].material;
            const auto mat_k = state.anode_particles[k].material;

            auto& mu_j = *state.anode_pairs.mu_pair_a[j][k];
            auto& mu_k = *state.anode_pairs.mu_pair_b[j][k];
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


static void InitializeAnodeParticles(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    const SimulationConfig& cfg, BoundaryConditions& bc, const std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_fields,
    const std::vector<double>& particle_totals, const std::vector<int>& particle_labels)
{
    const int np = static_cast<int>(particle_fields.size());
    state.anode_particles.clear();
    state.anode_particles.resize(np);

    // const std::vector<double>& init_values = cfg.init_anode_particles;

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
                    << " (label = " << particle_labels[k] << ")"
                    << std::endl;
        }

        auto& p = state.anode_particles[k];
        p.label = particle_labels[k];
        p.material = cfg.anode_materials[k];

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] Anode Particle " << k << " assigned material = ";

            switch (p.material)
            {
                case sim::MaterialType::Graphite:
                    std::cout << "Graphite";
                    break;

                case sim::MaterialType::Carbon:
                    std::cout << "Carbon";
                    break;
                
                default:
                {
                    mfem::mfem_error("Unsupported anode material chosen. Anode materials supported at this time: graphite, carbon.");
                }
            }

            std::cout << std::endl;
        }

        switch (p.material)
        {
            case sim::MaterialType::Graphite:
            {
                p.concentration = std::make_unique<ElectrodeCahnHilliard>(geometry, domain_parameters, p.material, cfg);
                break;
            }

            case sim::MaterialType::Carbon:
            {
                p.concentration = std::make_unique<ElectrodeDiffusion>(geometry, domain_parameters, p.material, cfg);
                break;
            }

            default:
            {
                mfem::mfem_error("Unsupported anode material physics.");
            }
        }

        p.Cn_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Cn_gf_psi     = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction      = std::make_unique<Reaction>(geometry, domain_parameters, cfg);
        p.Rxn_gf        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Rx_src        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        
        p.reaction->Initialize(*p.Rxn_gf, Constants::init_Rxn);

        const double init_cn = GetInitialValue(cfg.init_anode_particles, k, cfg.init_CnA);

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG]   Initial concentration = " << init_cn << std::endl;
        }

        p.concentration->SetupField(*p.Cn_gf, init_cn, *particle_fields[k], particle_totals[k]);
    
    }
}

static void InitializeCathodeParticles(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    const SimulationConfig& cfg, BoundaryConditions& bc, const std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_fields,
    const std::vector<double>& particle_totals, const std::vector<int>& particle_labels)
{
    const int np = static_cast<int>(particle_fields.size());
    state.cathode_particles.clear();

    if (np == 0)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout
                << "[DEBUG] No cathode particle groups were created from the geometry."
                << std::endl;
        }
        return;
    }

    if (cfg.cathode_materials.size() != static_cast<size_t>(np))
    {
        std::stringstream ss;
        ss << "Invalid number of cathode material entries.\n"
        << "The geometry contains " << np << " cathode particle group"
        << (np == 1 ? "" : "s")
        << ", but 'cathode_materials' contains "
        << cfg.cathode_materials.size() << " entr"
        << (cfg.cathode_materials.size() == 1 ? "y" : "ies") << ".\n\n"
        << "Please provide one material for each particle group.";
        mfem::mfem_error(ss.str().c_str());
    }

    if (cfg.init_cathode_particles.size() != static_cast<size_t>(np))
    {
        std::stringstream ss;
        ss << "Invalid number of cathode initial concentrations.\n"
        << "The geometry contains " << np << " cathode particle group"
        << (np == 1 ? "" : "s")
        << ", but 'init_cathode_particles' contains "
        << cfg.init_cathode_particles.size() << " entr"
        << (cfg.init_cathode_particles.size() == 1 ? "y" : "ies") << ".\n\n"
        << "Please provide one initial concentration for each particle group.";
        mfem::mfem_error(ss.str().c_str());
    }

    state.cathode_particles.resize(np);

    // const std::vector<double>& init_values = cfg.init_cathode_particles;

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
                    << " (label = " << particle_labels[k] << ")"
                    << std::endl;
        }

        auto& p = state.cathode_particles[k];
        p.label = particle_labels[k];
        p.material = cfg.cathode_materials[k];

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] Cathode Particle " << k << " assigned material = ";

            switch (p.material)
            {
                case sim::MaterialType::NMC:
                    std::cout << "NMC";
                    break;

                case sim::MaterialType::LFP:
                    std::cout << "LFP";
                    break;
                
                default:
                {
                    mfem::mfem_error("Unsupported cathode material chosen. Cathode materials supported at this time: NMC, LFP.");
                }
            }

            std::cout << std::endl;
        }

        switch (p.material)
        {
            case sim::MaterialType::NMC:
            {
                p.concentration = std::make_unique<ElectrodeDiffusion>(geometry, domain_parameters, p.material, cfg);
                break;
            }

            case sim::MaterialType::LFP:
            {
                p.concentration = std::make_unique<ElectrodeCahnHilliard>(geometry, domain_parameters, p.material, cfg);
                break;
            }

            default:
            {
                mfem::mfem_error("Unsupported cathode material physics.");
            }
        }

        p.Cn_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Cn_gf_psi     = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction      = std::make_unique<Reaction>(geometry, domain_parameters, cfg);
        p.Rxn_gf        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Rx_src        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction->Initialize(*p.Rxn_gf, Constants::init_Rxn);

        const double init_cn = GetInitialValue(cfg.init_cathode_particles, k, cfg.init_CnC);

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG]   Initial concentration = " << init_cn << std::endl;
        }

        p.concentration->SetupField(*p.Cn_gf, init_cn, *particle_fields[k], particle_totals[k]);
    }
}

void InitializeFields(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    BoundaryConditions& bc, const SimulationConfig& cfg)
{
    state.CnP_together = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.CnE_gf_psi   = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

    state.electrolyte_concentration = std::make_unique<ElectrolyteDiffusion>(geometry, domain_parameters, bc, cfg.mode, sim::MaterialType::Electrolyte, cfg);
    state.CnE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.electrolyte_concentration->SetupField(*state.CnE_gf, cfg.init_CnE, *domain_parameters.pse, domain_parameters.gtPse);
    
    state.electrolyte_potential = std::make_unique<ElectrolytePotential>(geometry, domain_parameters, bc, cfg.mode, cfg);
    state.phE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.electrolyte_potential->SetupField(*state.phE_gf, cfg.init_BvE, *domain_parameters.pse);

    state.reaction = std::make_unique<Reaction>(geometry, domain_parameters, cfg);
    state.Rxn_gf   = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.RxnA_gf  = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.RxnC_gf  = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.RxnE_gf  = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    
    state.reaction->Initialize(*state.Rxn_gf, Constants::init_Rxn);
    state.reaction->Initialize(*state.RxnA_gf, Constants::init_Rxn);
    state.reaction->Initialize(*state.RxnC_gf, Constants::init_Rxn);
    state.reaction->Initialize(*state.RxnE_gf, Constants::init_Rxn);

    if (cfg.mode == sim::CellMode::HALF)
    {
        if (cfg.half_electrode == sim::Electrode::ANODE)
        {
            const int np = static_cast<int>(domain_parameters.ps.size());

            state.CnA_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

            // state.anode_potential = std::make_unique<ElectrodePotential>(geometry, domain_parameters, bc, sim::Electrode::ANODE, sim::MaterialType::Graphite, cfg);
            MFEM_VERIFY(!cfg.anode_materials.empty(), "An anode material must be specified before initializing anode potential.");
            const sim::MaterialType anode_material = cfg.anode_materials.front();
            state.anode_potential = std::make_unique<ElectrodePotential>(geometry, domain_parameters, bc,sim::Electrode::ANODE, anode_material, cfg);

            state.phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.anode_potential->SetupField(*state.phA_gf, cfg.init_BvA, *domain_parameters.psi);

            InitializeAnodeParticles(state, geometry, domain_parameters, cfg, bc, domain_parameters.ps, domain_parameters.gtPs, domain_parameters.particle_labels);
            // InitializePairWorkspaces(state, geometry, static_cast<int>(state.anode_particles.size()));
            InitializePairWorkspaces(state.anode_pairs, geometry, static_cast<int>(state.anode_particles.size()),"anode");

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

            // state.cathode_potential = std::make_unique<ElectrodePotential>(geometry, domain_parameters, bc, sim::Electrode::CATHODE, sim::MaterialType::NMC, cfg);
            MFEM_VERIFY(!cfg.cathode_materials.empty(), "A cathode material must be specified before initializing cathode potential.");
            const sim::MaterialType cathode_material = cfg.cathode_materials.front();
            state.cathode_potential = std::make_unique<ElectrodePotential>(geometry, domain_parameters, bc, sim::Electrode::CATHODE, cathode_material, cfg);
            state.phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.cathode_potential->SetupField(*state.phC_gf, cfg.init_BvC, *domain_parameters.psi);

            InitializeCathodeParticles(state, geometry, domain_parameters, cfg, bc, domain_parameters.ps, domain_parameters.gtPs, domain_parameters.particle_labels);
            // InitializePairWorkspaces(state, geometry,static_cast<int>(state.cathode_particles.size()));
            InitializePairWorkspaces(state.cathode_pairs, geometry, static_cast<int>(state.cathode_particles.size()),"cathode");

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

        MFEM_VERIFY(!cfg.anode_materials.empty(), "An anode material must be specified before initializing anode potential.");
        const sim::MaterialType anode_material = cfg.anode_materials.front();
        MFEM_VERIFY(!cfg.cathode_materials.empty(), "A cathode material must be specified before initializing cathode potential.");
        const sim::MaterialType cathode_material = cfg.cathode_materials.front();

        state.anode_potential = std::make_unique<ElectrodePotential>(geometry, domain_parameters, bc, sim::Electrode::ANODE, anode_material, cfg);
        state.phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.anode_potential->SetupField(*state.phA_gf, cfg.init_BvA, *domain_parameters.psiA);

        state.cathode_potential = std::make_unique<ElectrodePotential>(geometry, domain_parameters, bc, sim::Electrode::CATHODE, cathode_material, cfg);
        state.phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.cathode_potential->SetupField(*state.phC_gf, cfg.init_BvC, *domain_parameters.psiC);

        InitializeAnodeParticles(state, geometry, domain_parameters, cfg, bc, domain_parameters.psA, domain_parameters.gtPsA, domain_parameters.anode_particle_labels);
        InitializeCathodeParticles(state, geometry, domain_parameters, cfg, bc, domain_parameters.psC, domain_parameters.gtPsC, domain_parameters.cathode_particle_labels);

        // InitializePairWorkspaces(state, geometry, static_cast<int>(state.anode_particles.size()));
        // InitializePairWorkspaces(state, geometry, static_cast<int>(state.cathode_particles.size()));

        const int npA = static_cast<int>(state.anode_particles.size());
        const int npC = static_cast<int>(state.cathode_particles.size());

        InitializePairWorkspaces(state.anode_pairs, geometry, npA, "anode");
        InitializePairWorkspaces(state.cathode_pairs, geometry, npC, "cathode");

        // Anode outputs
        state.anode_out.clear();
        state.anode_out.resize(npA);

        for (int k = 0; k < npA; ++k)
        {
            state.anode_out[k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        }

        // Cathode outputs
        state.cathode_out.clear();
        state.cathode_out.resize(npC);

        for (int k = 0; k < npC; ++k)
        {
            state.cathode_out[k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        }
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

void Pairs(PairWorkspaces& workspace, const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& weight_pairs,
    const std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>>& avp_pairs, int j,
    std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np, int t)
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

        MFEM_VERIFY(workspace.sum_pairs[a][b], "Pair sum workspace is null.");
        MFEM_VERIFY(workspace.mu_pair_a[a][b], "Pair chemical-potential workspace A is null.");
        MFEM_VERIFY(workspace.mu_pair_b[a][b], "Pair chemical-potential workspace B is null.");
        MFEM_VERIFY(weight_pairs[a][b], "Pair weight field is null.");
        MFEM_VERIFY(avp_pairs[a][b], "Pair interface field is null.");

        ConcentrationBase::PairCoupling pair;

        pair.sum_part = workspace.sum_pairs[a][b].get();
        pair.weight = weight_pairs[a][b].get();
        pair.grad_psi = avp_pairs[a][b].get();

        if (j < k)
        {
            pair.mu_self = workspace.mu_pair_a[a][b].get();
            pair.mu_nbr = workspace.mu_pair_b[a][b].get();
        }
        else
        {
            pair.mu_self = workspace.mu_pair_b[a][b].get();
            pair.mu_nbr = workspace.mu_pair_a[a][b].get();
        }

        pair_terms.push_back(pair);
    }
}
