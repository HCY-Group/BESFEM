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

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "\n===== Simulation Parameters =====\n"
                    << "output_dir = " << outdir << "\n"
                    << "dt   = " << cfg.dt   << "\n"
                    << "dh   = " << cfg.dh   << "\n"
                    << "gc   = " << cfg.gc   << "\n"
                    << "Cr   = " << cfg.Cr   << "\n"
                    << "Vsr0 = " << cfg.Vsr0 << "\n"
                    << "=================================\n"
                    << std::endl;
        }

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
            // std::cout << "Full cell mode is not yet implemented. Please use half-cell mode." << std::endl;
        }
        bc.SaveBoundaryConditionFields();

        // Define Adjuster for Surface Voltage & Current
        Adjust adjust(geometry, domain_parameters, cfg);

        // Initialize Concentration & Potential & Reaction Fields
        SimulationState state;
        InitializeFields(state, geometry, domain_parameters, bc, cfg);

        double VCell = 0.0;

        // ============================================================================
        // ===============================  TIME STEP LOOP  ===========================
        // ============================================================================

        if (cfg.mode == sim::CellMode::HALF)
        {
            
            // std::cout << "Starting simulation..." << std::endl;

            int t = 0;

            while (true) {

                double VCell = 0.0;

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    VCell = state.anode_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();
                }
                else
                {
                    VCell = state.cathode_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();
                }

                if (cfg.stop_mode == sim::StopMode::STEPS && t >= cfg.num_timesteps)
                {
                    break;
                }

                if (cfg.stop_mode == sim::StopMode::VOLTAGE && VCell <= cfg.VCut)
                {
                    break;
                }

            // for (int t = 0; t < cfg.num_timesteps; ++t) {

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    const int np = static_cast<int>(state.anode_particles.size());
                    std::vector<double> global_currents(np, 0.0);

                    UpdateAnodePairChemicalPotentials(state, geometry, domain_parameters.AvP_Pairs);

                    *state.Rxn_gf = 0.0;
                    for (int j = 0; j < np; ++j)
                    {
                        *state.anode_particles[j].Rx_src = *state.anode_particles[j].Rxn_gf;
                        *state.Rxn_gf += *state.anode_particles[j].Rxn_gf;

                        std::vector<ConcentrationBase::PairCoupling> pair_terms;
                        Pairs(state.anode_pairs, domain_parameters.WeightPairs, domain_parameters.AvP_Pairs, j, pair_terms, np, t);

                        state.anode_particles[j].concentration->UpdateConcentration(*state.anode_particles[j].Rx_src, *state.anode_particles[j].Cn_gf,
                            *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);

                    }

                    state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf,
                        *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});

                    if (t > 0 && t % 50 == 0) {
                        state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                    }

                    std::vector<mfem::ParGridFunction*> anode_cn_fields;
                    std::vector<mfem::ParGridFunction*> anode_ps_fields;
                    std::vector<sim::MaterialType> anode_materials;

                    anode_cn_fields.reserve(np);
                    anode_ps_fields.reserve(np);
                    anode_materials.reserve(np);

                    for (int j = 0; j < np; ++j)
                    {
                        anode_cn_fields.push_back(state.anode_particles[j].Cn_gf.get());
                        anode_ps_fields.push_back(domain_parameters.ps[j].get());
                        anode_materials.push_back(state.anode_particles[j].material);
                    }

                    state.anode_potential->AssembleSystem(anode_cn_fields, anode_ps_fields, anode_materials, *state.phA_gf);
                    state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    for (int j = 0; j < np; ++j)
                    {
                        state.anode_particles[j].reaction->ExchangeCurrentDensity(*state.anode_particles[j].Cn_gf, *domain_parameters.AvEs[j], state.anode_particles[j].material);
                    }

                    int iter = 0;
                    const int max_iter = 50;

                    double anode_time = 0.0;

                    while ((globalerror_P > 1e-6 || globalerror_E > 1e-6) && iter < max_iter) {
                        *state.Rxn_gf = 0.0;

                        for (int j = 0; j < np; ++j)
                        {
                            state.anode_particles[j].reaction->ButlerVolmer(*state.anode_particles[j].Rxn_gf, *state.anode_particles[j].Cn_gf,*state.CnE_gf,
                                *state.phA_gf, *state.phE_gf, *domain_parameters.AvEs[j]);
                            *state.Rxn_gf += *state.anode_particles[j].Rxn_gf;
                        }

                        state.anode_potential->UpdatePotential(*state.Rxn_gf, *state.phA_gf, *domain_parameters.psi, globalerror_P);
                        state.electrolyte_potential->UpdatePotential(*state.Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);

                        iter++;
                    }

                    if (iter == max_iter && mfem::Mpi::WorldRank() == 0) {
                        std::cout << "Warning: Maximum iterations reached at timestep " << t << " with Global Error P = " << globalerror_P << ", Global Error E = " << globalerror_E << std::endl;
                    }

                    for (int j = 0; j < np; ++j)
                    {
                        state.anode_particles[j].reaction->TotalReactionCurrent(*state.anode_particles[j].Rxn_gf, global_currents[j]);
                    }

                    double total_current = 0.0;
                    double total_target  = 0.0;

                    for (int j = 0; j < np; ++j)
                    {
                        total_current += global_currents[j];
                        total_target  += domain_parameters.gTrgPs[j];
                    }

                    double VCell = state.anode_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();

                    double sgn = std::copysign(1.0, total_target - total_current);
                    double dV  = cfg.dt * cfg.Vsr0 * sgn;

                    state.electrolyte_potential->AddBoundaryVoltage(dV);
                    *state.phE_gf += dV;


                    // ============================================================================
                    // ===============================  PRINT STATEMENTS  =========================
                    // ============================================================================

                    if (t % 5000 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::cout << "timestep: " << t << ", VCell = " << VCell << ", TotalCurrent = " << total_current << ", TotalTarget = " << total_target;

                        for (int j = 0; j < np; ++j)
                        {
                            std::cout << ", Current_" << j << " = " << global_currents[j] << ", Target_" << j << " = " << domain_parameters.gTrgPs[j];
                        }

                        std::cout << std::endl;
                    }

                    if (t % 5000 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        double XfrC_avg = 0.0;
                        double total_weight = 0.0;

                        std::cout << "timestep: " << t << " [ANODE HALF-CELL]" << ", VCell = " << VCell << ", BvE = " << state.electrolyte_potential->GetBoundaryVoltage();

                        std::cout
                        << "Cp_min = " << state.anode_particles[0].Cn_gf->Min()
                        << ", Cp_max = " << state.anode_particles[0].Cn_gf->Max()
                        << ", Ce_min = " << state.CnE_gf->Min()
                        << ", Ce_max = " << state.CnE_gf->Max()
                        << std::endl;

                        for (int j = 0; j < np; ++j)
                        {
                            const double Xfr_j = state.anode_particles[j].concentration->GetLithiation();
                            const double weight_j = domain_parameters.gtPs[j];

                            XfrC_avg += weight_j * Xfr_j;
                            total_weight += weight_j;

                            std::cout << ", Xfr_" << j << " = " << Xfr_j;
                        }

                        if (total_weight > 0.0)
                        {
                            XfrC_avg /= total_weight;
                        }

                        std::cout << ", XfrC_avg = " << XfrC_avg;
                        std::cout  << std::endl;
                    }
                    
                }
                    // ============================================================================
                    // ===============================  CATHODE  ==================================
                    // ============================================================================
                else
                {
                    // std::cout << "Updating cathode properties..." << std::endl;

                    const int np = static_cast<int>(state.cathode_particles.size());
                    std::vector<double> global_currents(np, 0.0);

                    UpdateCathodePairChemicalPotentials(state, geometry, domain_parameters.AvP_Pairs);

                    *state.Rxn_gf = 0.0;
                    for (int j = 0; j < np; ++j)
                    {
                        *state.cathode_particles[j].Rx_src = *state.cathode_particles[j].Rxn_gf;
                        *state.Rxn_gf += *state.cathode_particles[j].Rx_src;
                        
                        std::vector<ConcentrationBase::PairCoupling> pair_terms;
                        Pairs(state.cathode_pairs, domain_parameters.WeightPairs, domain_parameters.AvP_Pairs, j, pair_terms, np, t);

                        state.cathode_particles[j].concentration->UpdateConcentration(*state.cathode_particles[j].Rx_src, *state.cathode_particles[j].Cn_gf,
                            *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
                    }

                    // std::cout << "Cathode concentration updated." << std::endl;

                    // state.Rxn_gf->SaveAsOne("Rxn_after_concentration.gf");

                    state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf,
                        *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});
                        
                    if (t > 0 && t % 50 == 0) {
                        state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                    } 

                    // ============================================================
                    // Assemble one combined cathode potential
                    // ============================================================

                    if (t % 5 == 0)
                    {

                    std::vector<mfem::ParGridFunction*> cathode_cn_fields; // vector of pointers to cathode concentration fields
                    std::vector<mfem::ParGridFunction*> cathode_psi_fields; // vector of pointers to cathode potential fields
                    std::vector<sim::MaterialType> cathode_materials; // vector of cathode material types
 
                    cathode_cn_fields.reserve(np); // pre-allocate memory
                    cathode_psi_fields.reserve(np); // pre-allocate memory
                    cathode_materials.reserve(np); // pre-allocate memory

                    for (int j = 0; j < np; ++j)
                    {
                        cathode_cn_fields.push_back(state.cathode_particles[j].Cn_gf.get()); 
                        cathode_psi_fields.push_back(domain_parameters.ps[j].get());
                        cathode_materials.push_back(state.cathode_particles[j].material);
                    }

                    state.cathode_potential->AssembleSystem(cathode_cn_fields, cathode_psi_fields, cathode_materials, *state.phC_gf);
                    state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                    double globalerror_P = 1.0; // Error for particle potential
                    double globalerror_E = 1.0; // Error for electrolyte potential

                    for (int j = 0; j < np; ++j)
                    {
                        state.cathode_particles[j].reaction->ExchangeCurrentDensity(*state.cathode_particles[j].Cn_gf, *domain_parameters.AvEs[j], state.cathode_particles[j].material);
                    }

                        // while loop
                        int iter = 0;
                        const int max_iter = 50; // Maximum number of iterations to prevent infinite loops

                        // std::cout << "Starting iteration loop..." << std::endl;
                        while (globalerror_P > 1e-5 || globalerror_E > 1e-5 && iter < max_iter) {
                            *state.Rxn_gf = 0.0;

                            for (int j = 0; j < np; ++j) {
                                state.cathode_particles[j].reaction->ButlerVolmer(*state.cathode_particles[j].Rxn_gf, *state.cathode_particles[j].Cn_gf, *state.CnE_gf, *state.phC_gf, *state.phE_gf, *domain_parameters.AvEs[j]);
                                // CheckFinite(*state.cathode_particles[j].Rxn_gf, "particle Rxn");
                                *state.Rxn_gf += *state.cathode_particles[j].Rxn_gf;
                            }

                            // state.Rxn_gf->SaveAsOne("Rxn_TEST.gf");
                            // CheckFinite(*state.Rxn_gf, "total Rxn");

                            state.cathode_potential->UpdatePotential(*state.Rxn_gf, *state.phC_gf, *domain_parameters.psi, globalerror_P);
                            state.electrolyte_potential->UpdatePotential(*state.Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);

                            iter++;
                        }
                    }
                    
                    for (int j = 0; j < np; ++j){
                        state.cathode_particles[j].reaction->TotalReactionCurrent(*state.cathode_particles[j].Rxn_gf, global_currents[j]);
                    }

                    double total_current = 0.0;
                    double total_target = 0.0;

                    for (int j = 0; j < np; ++j)
                    {
                        total_current += global_currents[j];
                        total_target  += domain_parameters.gTrgPs[j];
                    }

                    double VCell = state.cathode_potential->GetBoundaryVoltage() - state.electrolyte_potential->GetBoundaryVoltage();

                    double sgn = std::copysign(1.0, total_target - total_current);
                    double dV  = cfg.dt * cfg.Vsr0 * sgn;

                    state.electrolyte_potential->AddBoundaryVoltage(dV);
                    *state.phE_gf += dV;


                    if (t % 5000 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::cout << "timestep: " << t << ", VCell = " << VCell << ", TotalCurrent = " << total_current << ", TotalTarget = " << total_target;

                        for (int j = 0; j < np; ++j)
                        {
                            std::cout << ", Current_" << j << " = " << global_currents[j] << ", Target_" << j << " = " << domain_parameters.gTrgPs[j];
                        }

                        std::cout << std::endl;
                    }

                    if (t % 5000 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        double XfrC_avg = 0.0;
                        double total_weight = 0.0;

                        std::cout << "timestep: " << t << " [CATHODE HALF-CELL]" << ", VCell = " << VCell << ", BvE = " << state.electrolyte_potential->GetBoundaryVoltage();

                        std::cout
                        << "Cp_min = " << state.cathode_particles[0].Cn_gf->Min()
                        << ", Cp_max = " << state.cathode_particles[0].Cn_gf->Max()
                        << ", Ce_min = " << state.CnE_gf->Min()
                        << ", Ce_max = " << state.CnE_gf->Max()
                        << std::endl;

                        for (int j = 0; j < np; ++j)
                        {
                            const double Xfr_j = state.cathode_particles[j].concentration->GetLithiation();
                            const double weight_j = domain_parameters.gtPs[j];

                            XfrC_avg += weight_j * Xfr_j;
                            total_weight += weight_j;

                            std::cout << ", Xfr_" << j << " = " << Xfr_j;
                        }

                        if (total_weight > 0.0)
                        {
                            XfrC_avg /= total_weight;
                        }

                        std::cout << ", XfrC_avg = " << XfrC_avg;
                        std::cout << std::endl;

                    }

                    // if (t % 10000 == 0)
                    // {
                    //     std::ostringstream rxn_path;
                    //     rxn_path << outdir << "/Rxn_" << t << ".gf";

                    //     state.Rxn_gf->SaveAsOne(rxn_path.str().c_str());
                    // }
                    
                }

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    std::vector<mfem::ParGridFunction*> anode_cn_fields;
                    anode_cn_fields.reserve(state.anode_particles.size());

                    for (auto &p : state.anode_particles)
                    {
                        anode_cn_fields.push_back(p.Cn_gf.get());
                    }

                    // Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters,
                    //     anode_cn_fields, state.anode_out, 10000);

                    Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters, anode_cn_fields, domain_parameters.ps, *domain_parameters.psi, state.anode_out, "A", 5000);
                }
                else
                {
                    std::vector<mfem::ParGridFunction*> cathode_cn_fields;
                    cathode_cn_fields.reserve(state.cathode_particles.size());

                    for (auto &p : state.cathode_particles)
                    {
                        cathode_cn_fields.push_back(p.Cn_gf.get());
                    }

                    // Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters,
                    //     cathode_cn_fields, state.cathode_out, 10000);
                    Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters, cathode_cn_fields, domain_parameters.ps, *domain_parameters.psi, state.cathode_out, "C", 5000);
                }

                t++;
            }

        } 
        else
        {
            // ============================================================================
            // ========================  FULL-CELL TIME STEPPING  ==========================
            // ============================================================================

            int t = 0;

            const int npA = static_cast<int>(state.anode_particles.size());
            const int npC = static_cast<int>(state.cathode_particles.size());

            MFEM_VERIFY(npA > 0, "Full-cell simulation contains no anode particles.");
            MFEM_VERIFY(npC > 0, "Full-cell simulation contains no cathode particles.");
            MFEM_VERIFY(state.RxnA_gf, "Full-cell anode reaction field RxnA_gf is not initialized.");
            MFEM_VERIFY(state.RxnC_gf, "Full-cell cathode reaction field RxnC_gf is not initialized.");

            if (mfem::Mpi::WorldRank() == 0)
            {
                std::cout
                    << "Starting full-cell simulation.\n"
                    << "    Anode particles:   " << npA << "\n"
                    << "    Cathode particles: " << npC
                    << std::endl;
            }

            while (true)
            {

                VCell = state.cathode_potential->GetBoundaryVoltage() - state.anode_potential->GetBoundaryVoltage();

                // ========================================================================
                // =========================  STOP CONDITIONS  =============================
                // ========================================================================

                if (cfg.stop_mode == sim::StopMode::STEPS && t >= cfg.num_timesteps)
                {
                    break;
                }

                if (cfg.stop_mode == sim::StopMode::VOLTAGE && VCell <= cfg.VCut)
                {
                    break;
                }

                // ========================================================================
                // ==============  UPDATE PARTICLE-PAIR CHEMICAL POTENTIALS  ===============
                // ========================================================================

                UpdateAnodePairChemicalPotentials(state, geometry, domain_parameters.AvP_PairsA);
                UpdateCathodePairChemicalPotentials(state, geometry, domain_parameters.AvP_PairsC);

                // ========================================================================
                // ====================  UPDATE ANODE CONCENTRATIONS  =======================
                // ========================================================================

                for (int j = 0; j < npA; ++j)
                {
                    MFEM_VERIFY(state.anode_particles[j].Rxn_gf, "Anode particle reaction field is null.");
                    MFEM_VERIFY(state.anode_particles[j].Rx_src, "Anode particle reaction source field is null.");
                    MFEM_VERIFY(state.anode_particles[j].Cn_gf, "Anode particle concentration field is null.");

                    // Freeze the previous reaction for the concentration solve.
                    *state.anode_particles[j].Rx_src = *state.anode_particles[j].Rxn_gf;

                    std::vector<ConcentrationBase::PairCoupling> pair_terms;

                    Pairs(state.anode_pairs, domain_parameters.WeightPairsA, domain_parameters.AvP_PairsA,
                        j, pair_terms, npA, t);

                    state.anode_particles[j].concentration->UpdateConcentration(*state.anode_particles[j].Rx_src,
                            *state.anode_particles[j].Cn_gf, *domain_parameters.psA[j], domain_parameters.gtPsA[j],
                            *domain_parameters.WeightEsA[j], pair_terms);
                }

                // ========================================================================
                // ===================  UPDATE CATHODE CONCENTRATIONS  ======================
                // ========================================================================

                for (int j = 0; j < npC; ++j)
                {
                    MFEM_VERIFY(state.cathode_particles[j].Rxn_gf, "Cathode particle reaction field is null.");
                    MFEM_VERIFY(state.cathode_particles[j].Rx_src, "Cathode particle reaction source field is null.");
                    MFEM_VERIFY(state.cathode_particles[j].Cn_gf, "Cathode particle concentration field is null.");

                    // Freeze the previous reaction for the concentration solve.
                    *state.cathode_particles[j].Rx_src = *state.cathode_particles[j].Rxn_gf;

                    std::vector<ConcentrationBase::PairCoupling> pair_terms;

                    Pairs(state.cathode_pairs, domain_parameters.WeightPairsC, domain_parameters.AvP_PairsC,
                        j, pair_terms, npC, t);

                    state.cathode_particles[j].concentration->UpdateConcentration(*state.cathode_particles[j].Rx_src,
                            *state.cathode_particles[j].Cn_gf, *domain_parameters.psC[j], domain_parameters.gtPsC[j],
                            *domain_parameters.WeightEsC[j], pair_terms);
                }

                // ========================================================================
                // ==============  BUILD ELECTRODE REACTION SOURCE FIELDS  ==================
                // ========================================================================

                *state.RxnA_gf = 0.0;
                *state.RxnC_gf = 0.0;

                for (int j = 0; j < npA; ++j)
                {
                    *state.RxnA_gf += *state.anode_particles[j].Rx_src;
                }

                for (int j = 0; j < npC; ++j)
                {
                    *state.RxnC_gf += *state.cathode_particles[j].Rx_src;
                }

                // // Optional combined field for output/debugging.
                // if (state.Rxn_gf)
                // {
                //     *state.Rxn_gf = *state.RxnA_gf;
                //     *state.Rxn_gf += *state.RxnC_gf;
                // }

                // ========================================================================
                // =================  UPDATE ELECTROLYTE CONCENTRATION  =====================
                // ========================================================================

                *state.RxnE_gf = 0.0;
                *state.RxnE_gf += *state.RxnA_gf;
                *state.RxnE_gf += *state.RxnC_gf;

                state.electrolyte_concentration->UpdateConcentration(*state.RxnE_gf, *state.CnE_gf, *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});

                if (t > 0 && t % 50 == 0)
                {
                    state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                }

                // ========================================================================
                // ====================  BUILD ANODE FIELD VECTORS  =========================
                // ========================================================================

                std::vector<mfem::ParGridFunction*>anode_cn_fields;
                std::vector<mfem::ParGridFunction*>anode_psi_fields;
                std::vector<sim::MaterialType>anode_materials;

                anode_cn_fields.reserve(npA);
                anode_psi_fields.reserve(npA);
                anode_materials.reserve(npA);

                for (int j = 0; j < npA; ++j)
                {
                    anode_cn_fields.push_back(state.anode_particles[j].Cn_gf.get());
                    anode_psi_fields.push_back(domain_parameters.psA[j].get());
                    anode_materials.push_back(state.anode_particles[j].material);
                }

                // ========================================================================
                // ===================  BUILD CATHODE FIELD VECTORS  ========================
                // ========================================================================

                std::vector<mfem::ParGridFunction*>cathode_cn_fields;
                std::vector<mfem::ParGridFunction*>cathode_psi_fields;
                std::vector<sim::MaterialType>cathode_materials;

                cathode_cn_fields.reserve(npC);
                cathode_psi_fields.reserve(npC);
                cathode_materials.reserve(npC);

                for (int j = 0; j < npC; ++j)
                {
                    cathode_cn_fields.push_back(state.cathode_particles[j].Cn_gf.get());
                    cathode_psi_fields.push_back(domain_parameters.psC[j].get());
                    cathode_materials.push_back(state.cathode_particles[j].material);
                }

                // ========================================================================
                // =====================  ASSEMBLE POTENTIAL SYSTEMS  =======================
                // ========================================================================

                state.anode_potential->AssembleSystem(anode_cn_fields, anode_psi_fields, anode_materials, *state.phA_gf);
                state.cathode_potential->AssembleSystem(cathode_cn_fields, cathode_psi_fields, cathode_materials, *state.phC_gf);
                state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

                // ========================================================================
                // =================  UPDATE EXCHANGE CURRENT DENSITIES  ====================
                // ========================================================================

                for (int j = 0; j < npA; ++j)
                {
                    state.anode_particles[j].reaction->ExchangeCurrentDensity(*state.anode_particles[j].Cn_gf,
                            *domain_parameters.AvEsA[j], state.anode_particles[j].material);
                }

                for (int j = 0; j < npC; ++j)
                {
                    state.cathode_particles[j].reaction->ExchangeCurrentDensity(*state.cathode_particles[j].Cn_gf,
                            *domain_parameters.AvEsC[j], state.cathode_particles[j].material);
                }

                // ========================================================================
                // ==============  COUPLED REACTION-POTENTIAL ITERATION  ====================
                // ========================================================================

                double globalerror_A = 1.0;
                double globalerror_C = 1.0;
                double globalerror_E = 1.0;

                int iter = 0;
                const int max_iter = 50;

                while ((globalerror_A > 1.0e-6 || globalerror_C > 1.0e-6 || globalerror_E > 1.0e-6) && iter < max_iter)
                {

                    *state.RxnA_gf = 0.0;
                    *state.RxnC_gf = 0.0;
                    *state.RxnE_gf = 0.0;

                    for (int j = 0; j < npA; ++j)
                    {
                        state.anode_particles[j].reaction->ButlerVolmer(*state.anode_particles[j].Rxn_gf,
                                *state.anode_particles[j].Cn_gf, *state.CnE_gf, *state.phA_gf, *state.phE_gf, *domain_parameters.AvEsA[j]);

                        *state.RxnA_gf += *state.anode_particles[j].Rxn_gf;
                    }

                    // --------------------------------------------------------------------
                    // Cathode Butler-Volmer reactions
                    // --------------------------------------------------------------------

                    for (int j = 0; j < npC; ++j)
                    {
                        state.cathode_particles[j].reaction->ButlerVolmer(*state.cathode_particles[j].Rxn_gf,
                                *state.cathode_particles[j].Cn_gf, *state.CnE_gf, *state.phC_gf, *state.phE_gf, *domain_parameters.AvEsC[j]);

                        *state.RxnC_gf += *state.cathode_particles[j].Rxn_gf;
                    }

                    // // Optional combined reaction for debugging/output.
                    // if (state.Rxn_gf)
                    // {
                    //     *state.Rxn_gf = *state.RxnA_gf;
                    //     *state.Rxn_gf += *state.RxnC_gf;
                    // }

                    *state.RxnE_gf = *state.RxnA_gf;
                    *state.RxnE_gf += *state.RxnC_gf;

                    // --------------------------------------------------------------------
                    // Update each solid potential using only its own reaction field.
                    // --------------------------------------------------------------------

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

                // ========================================================================
                // ======================  CALCULATE TOTAL CURRENTS  ========================
                // ========================================================================

                std::vector<double> anode_currents(npA, 0.0);
                std::vector<double> cathode_currents(npC, 0.0);

                for (int j = 0; j < npA; ++j)
                {
                    state.anode_particles[j].reaction->TotalReactionCurrent(*state.anode_particles[j].Rxn_gf, anode_currents[j]);
                }

                for (int j = 0; j < npC; ++j)
                {
                    state.cathode_particles[j].reaction->TotalReactionCurrent(*state.cathode_particles[j].Rxn_gf, cathode_currents[j]);
                }

                double global_current_A = 0.0;
                double global_current_C = 0.0;

                double total_target_A = 0.0;
                double total_target_C = 0.0;

                for (int j = 0; j < npA; ++j)
                {
                    global_current_A += anode_currents[j];
                    total_target_A += domain_parameters.gTrgPsA[j];
                }

                for (int j = 0; j < npC; ++j)
                {
                    global_current_C += cathode_currents[j];
                    total_target_C += domain_parameters.gTrgPsC[j];
                }

                // ========================================================================
                // ====================  ADJUST FULL-CELL CURRENT  ==========================
                // ========================================================================

                VCell = state.cathode_potential->GetBoundaryVoltage() - state.anode_potential->GetBoundaryVoltage();
                adjust.AdjustConstantCurrent(global_current_A, global_current_C, *state.anode_potential, *state.cathode_potential, *state.phA_gf, *state.phC_gf, VCell);
                VCell = state.cathode_potential->GetBoundaryVoltage() - state.anode_potential->GetBoundaryVoltage();

                // ========================================================================
                // =====================  CALCULATE LITHIATION  =============================
                // ========================================================================

                double XfrA_avg = 0.0;
                double XfrC_avg = 0.0;

                double total_anode_weight = 0.0;
                double total_cathode_weight = 0.0;

                for (int j = 0; j < npA; ++j)
                {
                    const double Xfr_j = state.anode_particles[j].concentration->GetLithiation();
                    const double weight_j = domain_parameters.gtPsA[j];
                    XfrA_avg += weight_j * Xfr_j;
                    total_anode_weight += weight_j;
                }

                for (int j = 0; j < npC; ++j)
                {
                    const double Xfr_j = state.cathode_particles[j].concentration->GetLithiation();
                    const double weight_j = domain_parameters.gtPsC[j];
                    XfrC_avg += weight_j * Xfr_j;
                    total_cathode_weight += weight_j;
                }

                if (total_anode_weight > 0.0)
                {
                    XfrA_avg /= total_anode_weight;
                }

                if (total_cathode_weight > 0.0)
                {
                    XfrC_avg /= total_cathode_weight;
                }

                // ========================================================================
                // ===========================  PRINT OUTPUT  ===============================
                // ========================================================================

                if (t % 1000 == 0 &&
                    mfem::Mpi::WorldRank() == 0)
                {
                    std::cout << "timestep: " << t
                        << " [FULL-CELL]"
                        << ", XfrA_avg = " << XfrA_avg
                        << ", XfrC_avg = " << XfrC_avg
                        << ", Anode current = " << global_current_A
                        << ", Cathode current = " << global_current_C
                        << ", Anode target = " << domain_parameters.gTrgI
                        << ", Cathode target = " << domain_parameters.gTrgI
                        << ", VCell = " << VCell
                        << ", BvA = " << state.anode_potential->GetBoundaryVoltage()
                        << ", BvC = " << state.cathode_potential->GetBoundaryVoltage()
                        << ", BvE = " << state.electrolyte_potential->GetBoundaryVoltage()
                        // << ", Ce_min = " << state.CnE_gf->Min()
                        // << ", Ce_max = " << state.CnE_gf->Max()
                        << std::endl;

                    for (int j = 0; j < npA; ++j)
                    {
                        std::cout << "    Anode particle " << j
                            // << ": current = " << anode_currents[j]
                            // << ", target = " << domain_parameters.gTrgPsA[j]
                            << ", Xfr = " << state.anode_particles[j].concentration->GetLithiation()
                            // << ", Cmin = " << state.anode_particles[j].Cn_gf->Min()
                            // << ", Cmax = " << state.anode_particles[j].Cn_gf->Max()
                            << std::endl;
                    }

                    for (int j = 0; j < npC; ++j)
                    {
                        std::cout << "    Cathode particle " << j
                            // << ": current = " << cathode_currents[j]
                            // << ", target = " << domain_parameters.gTrgPsC[j]
                            << ", Xfr = " << state.cathode_particles[j].concentration->GetLithiation()
                            // << ", Cmin = " << state.cathode_particles[j].Cn_gf->Min()
                            // << ", Cmax = " << state.cathode_particles[j].Cn_gf->Max()
                            << std::endl;
                    }
                }

                // ========================================================================
                // =============================  SAVE OUTPUT  ==============================
                // ========================================================================

                std::vector<mfem::ParGridFunction*>anode_cn_output;
                anode_cn_output.reserve(npA);

                for (auto& particle : state.anode_particles)
                {
                    anode_cn_output.push_back(particle.Cn_gf.get());
                }

                Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters, anode_cn_output, domain_parameters.psA, *domain_parameters.psiA, state.anode_out, "A", 1000);

                std::vector<mfem::ParGridFunction*>cathode_cn_output;
                cathode_cn_output.reserve(npC);

                for (auto& particle : state.cathode_particles)
                {
                    cathode_cn_output.push_back(particle.Cn_gf.get());
                }

                Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters, cathode_cn_output, domain_parameters.psC, *domain_parameters.psiC, state.cathode_out, "C", 1000);
                Utils::SaveCombinedElectrodeSnapshot(t, outdir, geometry, anode_cn_output, domain_parameters.psA, cathode_cn_output, domain_parameters.psC, 1000);

                // if (t % 1000 == 0)
                // {
                //     std::ostringstream rxnA_path;
                //     rxnA_path << outdir << "/RxnA_"  << t << ".gf";
                //     state.RxnA_gf->SaveAsOne(rxnA_path.str().c_str());

                //     std::ostringstream rxnC_path;
                //     rxnC_path << outdir << "/RxnC_" << t << ".gf"; 
                //     state.RxnC_gf->SaveAsOne(rxnC_path.str().c_str());
                // }

                ++t;
            }
        }

    }
    
    if (mfem::Mpi::WorldRank() == 0) { std::cout << "Simulation complete.\n"; }

    // End timing and output the total program execution time
    auto program_end = high_resolution_clock::now();
    if (mfem::Mpi::WorldRank() == 0) {std::cout << "Total Program Time: " 
            << duration_cast<seconds>(program_end - program_start).count() 
            << " seconds" << std::endl;}

    // Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    return 0;
}
