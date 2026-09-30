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
    Initialize_Geometry geometry(cfg, outdir);
    geometry.combine_particle_groups = cfg.combine_particle_groups;

    if (cfg.mode == sim::CellMode::HALF) {
        geometry.InitializeMesh(cfg.mesh_file, MPI_COMM_WORLD, cfg.order);
    } else {
        geometry.InitializeMesh(cfg.anode_mesh_file, cfg.cathode_mesh_file, MPI_COMM_WORLD, cfg.order);
    }

    // Initialize and Calculate Domain Parameters
    Domain_Parameters domain_parameters(geometry, cfg);
    domain_parameters.SetupDomainParameters(outdir);

    if (!cfg.geometry_only)
    {
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
        state.InitializeFields(geometry, domain_parameters, bc, cfg);

        double VCell = 0.0;

        // ============================================================================
        // =====================  HALF-CELL TIME STEP LOOP  ===========================
        // ============================================================================

        if (cfg.mode == sim::CellMode::HALF)
        {
            const bool is_anode = (cfg.half_electrode == sim::Electrode::ANODE);
            auto& electrode = is_anode ? state.anode : state.cathode;
            auto& particles = electrode.particles;
            auto& solid_potential = electrode.potential;
            auto& phS_gf = electrode.ph_gf;

            const int np = static_cast<int>(particles.size());

            double total_target = 0.0;
            for (int j = 0; j < np; ++j)
            {
                total_target += domain_parameters.gTrgPs[j];
            }

            int t = 0;

            double VCell_constant = 3.9;
            double VCell_difference = 0.0;

            while (true) {

                VCell = solid_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();
                // std::cout << "Timestep: " << t << ", VCell: " << VCell << std::endl;
                // std::cout << "Timestep: " << t << ", VCell_constant " << VCell_constant << std::endl;

                // VCell_difference = VCell - VCell_constant;
                // std::cout << "Timestep: " << t << ", VCell_difference " << VCell_difference << std::endl;

                if (Utils::ShouldStopSimulation(cfg, t, VCell)){break;}

                // PAIR CHEMICAL POTENTIALS
                electrode.UpdatePairChemicalPotentials(geometry, domain_parameters.AvP_Pairs);

                // PARTICLE CONCENTRATIONS
                electrode.UpdateParticleConcentrations(domain_parameters.WeightPairs, domain_parameters.AvP_Pairs, domain_parameters.ps, domain_parameters.gtPs, domain_parameters.WeightEs, *state.Rxn_gf);

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

                    electrode.BuildParticleFields(domain_parameters.ps, cn_fields, ps_fields, materials);

                    solid_potential->AssembleSystem(cn_fields, ps_fields, materials, *phS_gf);
                    state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                    electrode.UpdateExchangeCurrentDensity(domain_parameters.AvEs);

                    double globalerror_P = 1.0;
                    double globalerror_E = 1.0;

                    int iter = 0;
                    const int max_iter = 50;

                    while ((globalerror_P > 1e-5 || globalerror_E > 1e-5) && iter < max_iter)
                    {
                        electrode.UpdateButlerVolmerReactions(*state.Rxn_gf, *state.CnE_gf, *phS_gf, *state.phE_gf, domain_parameters.AvEs, domain_parameters.WeightEs);

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
                double total_current = electrode.CalculateElectrodeCurrent(global_currents);

                VCell = solid_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();

                if (cfg.Cr > 0 ? VCell >= VCell_constant : VCell <= VCell_constant)
                {
                    adjust.AdjustHalfCellCurrent(total_current, total_target, *state.electrolyte_potential, *state.phE_gf);
                }

                if (t % cfg.save_freq == 0)
                {
                    Utils::PrintHalfCellStatus(t, VCell, total_current, total_target, global_currents, state, domain_parameters, cfg.half_electrode);
                }

                Utils::SaveHalfCellSnapshot(t, outdir, geometry, domain_parameters, state, cfg.half_electrode, cfg.save_freq);

                ++t;
            }
        }
        // ============================================================================
        // ========================  FULL-CELL TIME STEPPING  =========================
        // ============================================================================
        else
        {
            int t = 0;

            const int npA = static_cast<int>(state.anode.particles.size());
            const int npC = static_cast<int>(state.cathode.particles.size());

            if (mfem::Mpi::WorldRank() == 0)
            {
                std::cout << "Starting full-cell simulation.\n" << "    Anode particles:   " << npA << "\n" << "    Cathode particles: " << npC << std::endl;
            }

            while (true)
            {

                VCell = state.cathode.potential->GetBoundaryVoltage() - state.anode.potential->GetBoundaryVoltage();

                if (Utils::ShouldStopSimulation(cfg, t, VCell)){break;}

                // PAIR CHEMICAL POTENTIALS
                state.anode.UpdatePairChemicalPotentials(geometry, domain_parameters.AvP_PairsA);
                state.cathode.UpdatePairChemicalPotentials(geometry, domain_parameters.AvP_PairsC);

                // PARTICLE CONCENTRATIONS
                state.anode.UpdateParticleConcentrations(domain_parameters.WeightPairsA, domain_parameters.AvP_PairsA, domain_parameters.psA, domain_parameters.gtPsA, domain_parameters.WeightEsA, *state.RxnA_gf);
                state.cathode.UpdateParticleConcentrations(domain_parameters.WeightPairsC, domain_parameters.AvP_PairsC, domain_parameters.psC, domain_parameters.gtPsC, domain_parameters.WeightEsC, *state.RxnC_gf);

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

                state.anode.BuildParticleFields(domain_parameters.psA, anode_cn_fields, anode_psi_fields, anode_materials);
                state.cathode.BuildParticleFields(domain_parameters.psC, cathode_cn_fields, cathode_psi_fields, cathode_materials);

                // ASSEMBLE POTENTIALS
                state.anode.potential->AssembleSystem(anode_cn_fields, anode_psi_fields, anode_materials, *state.anode.ph_gf);
                state.cathode.potential->AssembleSystem(cathode_cn_fields, cathode_psi_fields, cathode_materials, *state.cathode.ph_gf);
                state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                // EXCHANGE CURRENT DENSITY 
                state.anode.UpdateExchangeCurrentDensity(domain_parameters.AvEsA);
                state.cathode.UpdateExchangeCurrentDensity(domain_parameters.AvEsC);

                double globalerror_A = 1.0;
                double globalerror_C = 1.0;
                double globalerror_E = 1.0;

                int iter = 0;
                const int max_iter = 50;

                while ((globalerror_A > 1.0e-6 || globalerror_C > 1.0e-6 || globalerror_E > 1.0e-6) && iter < max_iter)
                {

                    state.anode.UpdateButlerVolmerReactions(*state.RxnA_gf, *state.CnE_gf, *state.anode.ph_gf, *state.phE_gf, domain_parameters.AvEsA, domain_parameters.WeightEsA);
                    state.cathode.UpdateButlerVolmerReactions(*state.RxnC_gf, *state.CnE_gf, *state.cathode.ph_gf, *state.phE_gf, domain_parameters.AvEsC, domain_parameters.WeightEsC);

                    *state.RxnE_gf = *state.RxnA_gf;
                    *state.RxnE_gf += *state.RxnC_gf;

                    state.anode.potential->UpdatePotential(*state.RxnA_gf, *state.anode.ph_gf, *domain_parameters.psiA, globalerror_A);
                    state.cathode.potential->UpdatePotential(*state.RxnC_gf, *state.cathode.ph_gf, *domain_parameters.psiC, globalerror_C);
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

                double global_current_A = state.anode.CalculateElectrodeCurrent(anode_currents);
                double global_current_C = state.cathode.CalculateElectrodeCurrent(cathode_currents);

                // ADJUST BOUNDARY VOLTAGES TO MAINTAIN GLOBAL CURRENT CONSERVATION
                VCell = state.cathode.potential->GetBoundaryVoltage() - state.anode.potential->GetBoundaryVoltage();
                adjust.AdjustConstantCurrent(global_current_A, global_current_C, *state.anode.potential, *state.cathode.potential, *state.anode.ph_gf, *state.cathode.ph_gf, VCell);
                VCell = state.cathode.potential->GetBoundaryVoltage() - state.anode.potential->GetBoundaryVoltage();

                if (t % cfg.save_freq == 0)
                {
                    Utils::PrintFullCellStatus(t, VCell, global_current_A, global_current_C, state, domain_parameters);
                }

                Utils::SaveFullCellSnapshot(t, outdir, geometry, domain_parameters, state, cfg.save_freq);

                ++t;
            }
        }
    }

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << (cfg.geometry_only ? "Geometry-only setup complete. Fields saved to "
                                       : "Simulation complete. Output saved to ")
                  << outdir << '\n';
    }
    }
    

    auto program_end = std::chrono::high_resolution_clock::now();

    Utils::PrintProgramTime(program_start, program_end);

    mfem::Hypre::Finalize();
    mfem::Mpi::Finalize();

    return 0;
}
