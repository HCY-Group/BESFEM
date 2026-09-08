#include "mfem.hpp"
#include "mpi.h"
#include "../include/BESFEM_All.hpp"

#include <chrono>
#include <iostream>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <vector>

int main(int argc, char *argv[]) {

    // Start measuring the program execution time
    using namespace std::chrono;
    auto program_start = high_resolution_clock::now();

    // Initialize MPI for parallel processing and HYPRE for solver setup
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    {

    SimulationConfig cfg = ParseSimulationArgs(argc, argv);
    ValidateConfig(cfg, argc, argv);

    std::string outdir = Utils::BuildRunOutdir();
    if (mfem::Mpi::WorldRank() == 0)
    {
        std::filesystem::create_directories(outdir);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // ============================================================================
    // ===============================  START SIMULATION  =========================
    // ============================================================================

    Utils::PrintSimulationParameters(cfg, outdir);

    // Initialize Mesh & Geometry
    Initialize_Geometry geometry(cfg);
    geometry.combine_particle_groups = cfg.combine_particle_groups;

    if (cfg.mode == sim::CellMode::HALF) {
        geometry.InitializeMesh(cfg.mesh_file, MPI_COMM_WORLD, cfg.order);
    } else {
        geometry.InitializeMesh(cfg.anode_mesh_file, cfg.cathode_mesh_file, MPI_COMM_WORLD, cfg.order);
    }

    // Initialize and Calculate Domain Parameters
    Domain_Parameters domain_parameters(geometry, cfg);
    domain_parameters.SetupDomainParameters();

    // Initialize Boundary Conditions 
    BoundaryConditions bc(geometry, domain_parameters);
    if (cfg.mode == sim::CellMode::HALF) {
        bc.SetupBoundaryConditions(sim::CellMode::HALF, cfg.half_electrode);
    } else {
        bc.SetupBoundaryConditions(sim::CellMode::FULL, sim::Electrode::BOTH);
    }
    bc.SaveBoundaryConditionFields();

    // Define Adjuster for Surface Voltage & Current
    Adjust adjust(geometry, domain_parameters, cfg);

    // Initialize Concentration & Potential & Reaction Fields
    SimulationState state;
    InitializeFields(state, geometry, domain_parameters, bc, cfg);

    double VCell = 0.0;

    // ============================================================================
    // =====================  HALF-CELL TIME STEP LOOP  ===========================
    // ============================================================================

        if (cfg.mode == sim::CellMode::HALF)
        {
            const bool is_anode = (cfg.half_electrode == sim::Electrode::ANODE);
            auto& particles = is_anode ? state.anode_particles : state.cathode_particles;
            auto& pairs = is_anode ? state.anode_pairs : state.cathode_pairs;
            auto& solid_potential = is_anode ? state.anode_potential : state.cathode_potential;
            auto& phS_gf = is_anode ? state.phA_gf : state.phC_gf;

            const int np = static_cast<int>(particles.size());

            double total_target = 0.0;
            for (int j = 0; j < np; ++j)
            {
                total_target += domain_parameters.gTrgPs[j];
            }

            int t = 0;
    
            while (true) {

                VCell = solid_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();

                if (Utils::ShouldStopSimulation(cfg, t, VCell)){break;}

                // PAIR CHEMICAL POTENTIALS
                UpdatePairChemicalPotentials(particles, pairs, geometry, domain_parameters.AvP_Pairs);

                // PARTICLE CONCENTRATIONS
                UpdateParticleConcentrations(particles, pairs, domain_parameters.WeightPairs, domain_parameters.AvP_Pairs, domain_parameters.ps, domain_parameters.gtPs, domain_parameters.WeightEs, *state.Rxn_gf, t);

                // ELECTROLYTE CONCENTRATION
                state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf, *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});

                if (t > 0 && t % 50 == 0)
                {
                    state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                }

                // POTENTIALS
                if (t % 5 == 0)
                {
                    std::vector<mfem::ParGridFunction*> cn_fields;
                    std::vector<mfem::ParGridFunction*> ps_fields;
                    std::vector<sim::MaterialType> materials;

                    BuildParticleFields(particles, domain_parameters.ps, cn_fields, ps_fields, materials);

                    solid_potential->AssembleSystem(cn_fields, ps_fields, materials, *phS_gf);
                    state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                    UpdateExchangeCurrentDensity(particles, domain_parameters.AvEs);

                    double globalerror_P = 1.0;
                    double globalerror_E = 1.0;

                    int iter = 0;
                    const int max_iter = 50;

                    while ((globalerror_P > 1e-5 || globalerror_E > 1e-5) && iter < max_iter)
                    {
                        UpdateButlerVolmerReactions(particles, *state.Rxn_gf, *state.CnE_gf, *phS_gf, *state.phE_gf, domain_parameters.AvEs, domain_parameters.WeightEs);

                        solid_potential->UpdatePotential(*state.Rxn_gf, *phS_gf, *domain_parameters.psi, globalerror_P);
                        state.electrolyte_potential->UpdatePotential(*state.Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);

                        ++iter;
                    }

                    if (iter == max_iter && mfem::Mpi::WorldRank() == 0)
                    { 
                        std::cout << "Warning: half-cell potential iteration reached " << max_iter << " iterations at timestep " << t << ". Error_P = " << globalerror_P << ", Error_E = " << globalerror_E << std::endl;
                    }
                }

                std::vector<double> global_currents;
                double total_current = CalculateElectrodeCurrent(particles, global_currents);

                VCell = solid_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();

                adjust.AdjustHalfCellCurrent(total_current, total_target, *state.electrolyte_potential, *state.phE_gf);

                if (t % 5000 == 0)
                {
                    Utils::PrintHalfCellStatus(t, VCell, total_current, total_target, global_currents, state, domain_parameters, cfg.half_electrode);
                }

                Utils::SaveHalfCellSnapshot(t, outdir, geometry, domain_parameters, state, cfg.half_electrode, 5000);

                ++t;
            }
        }
        // ============================================================================
        // ========================  FULL-CELL TIME STEPPING  =========================
        // ============================================================================
        else
        {
            int t = 0;

            const int npA = static_cast<int>(state.anode_particles.size());
            const int npC = static_cast<int>(state.cathode_particles.size());

            if (mfem::Mpi::WorldRank() == 0)
            {
                std::cout << "Starting full-cell simulation.\n" << "    Anode particles:   " << npA << "\n" << "    Cathode particles: " << npC << std::endl;
            }

            while (true)
            {

                VCell = state.cathode_potential->GetBoundaryVoltage() - state.anode_potential->GetBoundaryVoltage();

                if (Utils::ShouldStopSimulation(cfg, t, VCell)){break;}

                // PAIR CHEMICAL POTENTIALS
                UpdatePairChemicalPotentials(state.anode_particles, state.anode_pairs, geometry, domain_parameters.AvP_PairsA);
                UpdatePairChemicalPotentials(state.cathode_particles, state.cathode_pairs, geometry, domain_parameters.AvP_PairsC);

                // PARTICLE CONCENTRATIONS
                UpdateParticleConcentrations(state.anode_particles, state.anode_pairs, domain_parameters.WeightPairsA, domain_parameters.AvP_PairsA, domain_parameters.psA, domain_parameters.gtPsA, domain_parameters.WeightEsA, *state.RxnA_gf, t);
                UpdateParticleConcentrations(state.cathode_particles, state.cathode_pairs, domain_parameters.WeightPairsC, domain_parameters.AvP_PairsC, domain_parameters.psC, domain_parameters.gtPsC, domain_parameters.WeightEsC, *state.RxnC_gf, t);

                *state.RxnE_gf = 0.0;
                *state.RxnE_gf += *state.RxnA_gf;
                *state.RxnE_gf += *state.RxnC_gf;

                // ELECTROLYTE CONCENTRATION
                state.electrolyte_concentration->UpdateConcentration(*state.RxnE_gf, *state.CnE_gf, *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});

                if (t > 0 && t % 50 == 0)
                {
                    state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                }

                std::vector<mfem::ParGridFunction*> anode_cn_fields;
                std::vector<mfem::ParGridFunction*> anode_psi_fields;
                std::vector<sim::MaterialType> anode_materials;

                std::vector<mfem::ParGridFunction*> cathode_cn_fields;
                std::vector<mfem::ParGridFunction*> cathode_psi_fields;
                std::vector<sim::MaterialType> cathode_materials;

                BuildParticleFields(state.anode_particles, domain_parameters.psA, anode_cn_fields, anode_psi_fields, anode_materials);
                BuildParticleFields(state.cathode_particles, domain_parameters.psC, cathode_cn_fields, cathode_psi_fields, cathode_materials);

                // ASSEMBLE POTENTIALS
                state.anode_potential->AssembleSystem(anode_cn_fields, anode_psi_fields, anode_materials, *state.phA_gf);
                state.cathode_potential->AssembleSystem(cathode_cn_fields, cathode_psi_fields, cathode_materials, *state.phC_gf);
                state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                // EXCHANGE CURRENT DENSITY 
                UpdateExchangeCurrentDensity(state.anode_particles, domain_parameters.AvEsA);
                UpdateExchangeCurrentDensity(state.cathode_particles, domain_parameters.AvEsC);

                double globalerror_A = 1.0;
                double globalerror_C = 1.0;
                double globalerror_E = 1.0;

                int iter = 0;
                const int max_iter = 50;

                while ((globalerror_A > 1.0e-6 || globalerror_C > 1.0e-6 || globalerror_E > 1.0e-6) && iter < max_iter)
                {

                    UpdateButlerVolmerReactions(state.anode_particles, *state.RxnA_gf, *state.CnE_gf, *state.phA_gf, *state.phE_gf, domain_parameters.AvEsA, domain_parameters.WeightEsA);
                    UpdateButlerVolmerReactions(state.cathode_particles, *state.RxnC_gf, *state.CnE_gf, *state.phC_gf, *state.phE_gf, domain_parameters.AvEsC, domain_parameters.WeightEsC);

                    *state.RxnE_gf = *state.RxnA_gf;
                    *state.RxnE_gf += *state.RxnC_gf;

                    state.anode_potential->UpdatePotential(*state.RxnA_gf, *state.phA_gf, *domain_parameters.psiA, globalerror_A);
                    state.cathode_potential->UpdatePotential(*state.RxnC_gf, *state.phC_gf, *domain_parameters.psiC, globalerror_C);
                    state.electrolyte_potential->UpdatePotential(*state.RxnE_gf, *state.phE_gf, *domain_parameters.pse,  globalerror_E);

                    ++iter;
                }

                if (iter == max_iter && mfem::Mpi::WorldRank() == 0)
                {
                    std::cout << "Warning: full-cell potential iteration reached " << max_iter << " iterations at timestep " << t
                        << ". Error_A = " << globalerror_A << ", Error_C = " << globalerror_C << ", Error_E = " << globalerror_E
                        << std::endl;
                }

                std::vector<double> anode_currents;
                std::vector<double> cathode_currents;

                double global_current_A = CalculateElectrodeCurrent(state.anode_particles, anode_currents);
                double global_current_C = CalculateElectrodeCurrent(state.cathode_particles, cathode_currents);

                // ADJUST BOUNDARY VOLTAGES TO MAINTAIN GLOBAL CURRENT CONSERVATION
                VCell = state.cathode_potential->GetBoundaryVoltage() - state.anode_potential->GetBoundaryVoltage();
                adjust.AdjustConstantCurrent(global_current_A, global_current_C, *state.anode_potential, *state.cathode_potential, *state.phA_gf, *state.phC_gf, VCell);
                VCell = state.cathode_potential->GetBoundaryVoltage() - state.anode_potential->GetBoundaryVoltage();

                if (t % 1000 == 0)
                {
                    Utils::PrintFullCellStatus(t, VCell, global_current_A, global_current_C, state, domain_parameters);
                }

                Utils::SaveFullCellSnapshot(t, outdir, geometry, domain_parameters, state, 1000);

                ++t;
            }
        }
    }
    
    if (mfem::Mpi::WorldRank() == 0) { std::cout << "Simulation complete.\n"; }

    auto program_end = std::chrono::high_resolution_clock::now();

    Utils::PrintProgramTime(program_start, program_end);

    mfem::Hypre::Finalize();
    mfem::Mpi::Finalize();

    return 0;
}
