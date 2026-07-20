// Unit tests and benchmark tests
// compile with make test
// run with: ./unit_tests -cfg ../tests/test_run_config.txt

// Benchmark tests implemented:
// -ElectrolytePotential::UpdatePotential
// -ElectrodePotential::UpdatePotential

// Benchmark tests planned (but not yet implemented):
// -ElectrolyteDiffusion::UpdateConcentration
// -ElectrodeDiffusion::UpdateConcentration
// -ElectrodeCahnHilliard::UpdateConcentration

// Other Tests we should do, need to plan how to do them:
// -TiffReader
// -Test if "islands" are still in geometry
// -probably at least one for each source file...



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

    MPI_Barrier(MPI_COMM_WORLD);


        // Initialize Mesh & Geometry
        Initialize_Geometry geometry(cfg);
        geometry.combine_particle_groups = cfg.combine_particle_groups;

        if (cfg.mode == sim::CellMode::HALF) {
            geometry.InitializeMesh(cfg.mesh_file, MPI_COMM_WORLD, cfg.order);
        } else {
            "Full cell mode is not yet implemented. Please use half-cell mode.";
        }

        // Initialize and Calculate Domain Parameters
        Domain_Parameters domain_parameters(geometry, cfg);
        domain_parameters.SetupDomainParameters();

        // Initialize Boundary Conditions 
        BoundaryConditions bc(geometry, domain_parameters);
        if (cfg.mode == sim::CellMode::HALF) {
            bc.SetupBoundaryConditions(sim::CellMode::HALF, cfg.half_electrode);
        } else {
            // bc.SetupBoundaryConditions(sim::CellMode::FULL, sim::Electrode::BOTH);
            "Full cell mode is not yet implemented. Please use half-cell mode.";
        }

        // Define Adjuster for Surface Voltage & Current
        Adjust adjust(geometry, domain_parameters, cfg);

        // Initialize Concentration & Potential & Reaction Fields
        SimulationState state;
        InitializeFields(state, geometry, domain_parameters, bc, cfg);



        // =====================
        // poisson equation
        // =====================
        if (cfg.half_electrode == sim::Electrode::ANODE)
        {
        std::vector<mfem::ParGridFunction*> anode_cn_fields; // vector of pointers to cathode concentration fields
        std::vector<mfem::ParGridFunction*> anode_psi_fields; // vector of pointers to cathode potential fields
        std::vector<sim::MaterialType> anode_materials; // vector of cathode material types
                    
        const int np = static_cast<int>(state.anode_particles.size());
 
        anode_cn_fields.reserve(np); // pre-allocate memory
        anode_psi_fields.reserve(np); // pre-allocate memory
        anode_materials.reserve(np); // pre-allocate memory

        for (int j = 0; j < np; ++j)
        {
            anode_cn_fields.push_back(state.anode_particles[j].Cn_gf.get()); 
            anode_psi_fields.push_back(domain_parameters.ps[j].get());
            anode_materials.push_back(state.anode_particles[j].material);
        }
        state.anode_potential->AssembleSystem(anode_cn_fields, anode_psi_fields, anode_materials, *state.phA_gf);
        state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

        double globalerror_P = 1.0; // Error for particle potential
        double globalerror_E = 1.0; // Error for electrolyte potential

        mfem::ParGridFunction Rxn_p(*domain_parameters.AvEs[0]);
        mfem::ParGridFunction Rxn_e(*domain_parameters.AvEs[0]);

        //state.cathode_potential->UpdatePotential(*state.Rxn_gf, *state.phC_gf, *domain_parameters.psi, globalerror_P);
        //state.electrolyte_potential->UpdatePotential(*state.Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);
        state.anode_potential->UpdatePotential(Rxn_p, *state.phA_gf, *domain_parameters.psi, globalerror_P);
        state.electrolyte_potential->UpdatePotential(Rxn_e, *state.phE_gf, *domain_parameters.psi, globalerror_E);
        
        state.Rxn_gf->SaveAsOne("rxn_test");
        domain_parameters.psi->SaveAsOne("psi_test");
        domain_parameters.pse->SaveAsOne("pse_test");
        state.phA_gf->SaveAsOne("phiA");
        state.phE_gf->SaveAsOne("phiE");
        std::cout << "Faraday: " << Constants::Frd << std::endl;
        //TODO: figure out how to return kappa values from electrode and electrolyte 
        //std::cout << "Global errors: " << globalerror_P << " " << globalerror_E << std::endl;

        // ANALYTICAL SOLUTION
        Rxn_p.FESpace()->GetMesh()->EnsureNodes();
        mfem::GridFunction *nodes = Rxn_p.FESpace()->GetMesh()->GetNodes();
        //std::cout << "node 0: " << *nodes[0] << std::endl; 
        //std::cout << "node 1: " << *nodes[1] << std::endl;
        mfem::GridFunction x;
        mfem::GridFunction y;
        x.MakeRef(*nodes, /*offset=*/0,                                   /*size=*/Rxn_p.FESpace()->GetMesh()->GetNV());
        y.MakeRef(*nodes, /*offset=*/Rxn_p.FESpace()->GetMesh()->GetNV(), /*size=*/Rxn_p.FESpace()->GetMesh()->GetNV());
        //std::cout << "x min and max: " << x.Min() << " " << x.Max() << std::endl; 
        //std::cout << "y min and max: " << y.Min() << " " << y.Max() << std::endl; 
        
        mfem::ParGridFunction phA_an(*state.phA_gf);
        std::cout << "kappa: " << state.anode_potential->GetConductivity().Min() << " " << 
                                  state.anode_potential->GetConductivity().Max() << std::endl;
        std::cout << "F/kap: " << Constants::Frd/state.anode_potential->GetConductivity().Max() << std::endl;

        /*for (i = 0; i < phA_an.Size(); i++) 
        {
            phA_an(i) = y(i)*
        }*/

        } 
        
        else
        {
        std::vector<mfem::ParGridFunction*> cathode_cn_fields; // vector of pointers to cathode concentration fields
        std::vector<mfem::ParGridFunction*> cathode_psi_fields; // vector of pointers to cathode potential fields
        std::vector<sim::MaterialType> cathode_materials; // vector of cathode material types
                    
        const int np = static_cast<int>(state.cathode_particles.size());
 
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

        mfem::ParGridFunction Rxn_p(*domain_parameters.AvEs[0]);
        mfem::ParGridFunction Rxn_e(*domain_parameters.AvEs[0]);

        //state.cathode_potential->UpdatePotential(*state.Rxn_gf, *state.phC_gf, *domain_parameters.psi, globalerror_P);
        //state.electrolyte_potential->UpdatePotential(*state.Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);
        state.cathode_potential->UpdatePotential(Rxn_p, *state.phC_gf, *domain_parameters.psi, globalerror_P);
        state.electrolyte_potential->UpdatePotential(Rxn_e, *state.phE_gf, *domain_parameters.psi, globalerror_E);
        
        state.Rxn_gf->SaveAsOne("rxn_test");
        domain_parameters.psi->SaveAsOne("psi_test");
        domain_parameters.pse->SaveAsOne("pse_test");
        state.phC_gf->SaveAsOne("phiC");
        state.phE_gf->SaveAsOne("phiE");
        //std::cout << "Global errors: " << globalerror_P << " " << globalerror_E << std::endl;

        // ANALYTICAL SOLUTION
        Rxn_p.FESpace()->GetMesh()->EnsureNodes();
        mfem::GridFunction *nodes = Rxn_p.FESpace()->GetMesh()->GetNodes();
        //std::cout << "node 0: " << *nodes[0] << std::endl; 
        //std::cout << "node 1: " << *nodes[1] << std::endl;
        mfem::GridFunction x;
        mfem::GridFunction y;
        x.MakeRef(*nodes, /*offset=*/0,                                   /*size=*/Rxn_p.FESpace()->GetMesh()->GetNV());
        y.MakeRef(*nodes, /*offset=*/Rxn_p.FESpace()->GetMesh()->GetNV(), /*size=*/Rxn_p.FESpace()->GetMesh()->GetNV());
        //std::cout << "x min and max: " << x.Min() << " " << x.Max() << std::endl; 
        //std::cout << "y min and max: " << y.Min() << " " << y.Max() << std::endl; 
        
        mfem::ParGridFunction phC_an(*state.phC_gf);
        std::cout << "kappa: " << state.cathode_potential->GetConductivity().Min() << " " << 
                                  state.cathode_potential->GetConductivity().Max() << std::endl;







        }
/*
        // double VCell = 0.0;

        // ============================================================================
        // ===============================  TIME STEP LOOP  ===========================
        // ============================================================================

        if (cfg.mode == sim::CellMode::HALF)
        {
            
            int t = 0;

            while (true) {

                double VCell = 0.0;

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    VCell = state.anode_potential->GetBoundaryVoltage()
                        - state.electrolyte_potential->GetBoundaryVoltage();
                }
                else
                {
                    VCell = state.cathode_potential->GetBoundaryVoltage()
                        - state.electrolyte_potential->GetBoundaryVoltage();
                }

                if (cfg.stop_mode == sim::StopMode::STEPS &&
                    t >= cfg.num_timesteps)
                {
                    break;
                }

                if (cfg.stop_mode == sim::StopMode::VOLTAGE &&
                    VCell <= cfg.VCut)
                {
                    break;
                }

            // for (int t = 0; t < cfg.num_timesteps; ++t) {

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    const int np = static_cast<int>(state.anode_particles.size());
                    std::vector<double> global_currents(np, 0.0);

                    UpdateAnodePairChemicalPotentials(state, geometry, domain_parameters);

                    *state.Rxn_gf = 0.0;
                    for (int j = 0; j < np; ++j)
                    {
                        *state.anode_particles[j].Rx_src = *state.anode_particles[j].Rxn_gf;
                        *state.Rxn_gf += *state.anode_particles[j].Rxn_gf;

                        std::vector<ConcentrationBase::PairCoupling> pair_terms;
                        Pairs(state, geometry, domain_parameters, j, pair_terms, np, t);

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

                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::cout << "timestep: " << t << ", VCell = " << VCell << ", TotalCurrent = " << total_current << ", TotalTarget = " << total_target;

                        for (int j = 0; j < np; ++j)
                        {
                            std::cout << ", Current_" << j << " = " << global_currents[j] << ", Target_" << j << " = " << domain_parameters.gTrgPs[j];
                        }

                        std::cout << std::endl;
                    }

                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
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
                    const int np = static_cast<int>(state.cathode_particles.size());
                    std::vector<double> global_currents(np, 0.0);

                    UpdateCathodePairChemicalPotentials(state, geometry, domain_parameters);

                    *state.Rxn_gf = 0.0;
                    for (int j = 0; j < np; ++j)
                    {
                        *state.cathode_particles[j].Rx_src = *state.cathode_particles[j].Rxn_gf;
                        *state.Rxn_gf += *state.cathode_particles[j].Rxn_gf;
                        
                        std::vector<ConcentrationBase::PairCoupling> pair_terms;
                        Pairs(state, geometry, domain_parameters, j, pair_terms, np, t);
                        
                        state.cathode_particles[j].concentration->UpdateConcentration(*state.cathode_particles[j].Rx_src, *state.cathode_particles[j].Cn_gf,
                            *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
                    }

                    state.Rxn_gf->SaveAsOne("Rxn_after_concentration.gf");

                    state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf,
                        *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});
                        
                    if (t > 0 && t % 50 == 0) {
                        state.electrolyte_concentration->SaltConservation(*state.CnE_gf, *domain_parameters.pse);
                    } 

                    // ============================================================
                    // Assemble one combined cathode potential
                    // ============================================================

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

                        while (globalerror_P > 1e-5 || globalerror_E > 1e-5 && iter < max_iter) {
                            *state.Rxn_gf = 0.0;

                            for (int j = 0; j < np; ++j) {
                                state.cathode_particles[j].reaction->ButlerVolmer(*state.cathode_particles[j].Rxn_gf, *state.cathode_particles[j].Cn_gf, *state.CnE_gf, *state.phC_gf, *state.phE_gf, *domain_parameters.AvEs[j]);
                                *state.Rxn_gf += *state.cathode_particles[j].Rxn_gf;
                            }
                            state.cathode_potential->UpdatePotential(*state.Rxn_gf, *state.phC_gf, *domain_parameters.psi, globalerror_P);
                            state.electrolyte_potential->UpdatePotential(*state.Rxn_gf, *state.phE_gf, *domain_parameters.pse, globalerror_E);

                            iter++;

                        }

                        if (iter == max_iter && mfem::Mpi::WorldRank() == 0) {
                            std::cout << "Warning: Maximum iterations reached at timestep " << t << " with Global Error P = " << globalerror_P << ", Global Error E = " << globalerror_E << std::endl;
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


                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
                    {
                        std::cout << "timestep: " << t << ", VCell = " << VCell << ", TotalCurrent = " << total_current << ", TotalTarget = " << total_target;

                        for (int j = 0; j < np; ++j)
                        {
                            std::cout << ", Current_" << j << " = " << global_currents[j] << ", Target_" << j << " = " << domain_parameters.gTrgPs[j];
                        }

                        std::cout << std::endl;
                    }

                    if (t % 100 == 0 && mfem::Mpi::WorldRank() == 0)
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
                    
                }

                if (cfg.half_electrode == sim::Electrode::ANODE)
                {
                    std::vector<mfem::ParGridFunction*> anode_cn_fields;
                    anode_cn_fields.reserve(state.anode_particles.size());

                    for (auto &p : state.anode_particles)
                    {
                        anode_cn_fields.push_back(p.Cn_gf.get());
                    }

                    Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters,
                        anode_cn_fields, state.anode_out, 100);
                }
                else
                {
                    std::vector<mfem::ParGridFunction*> cathode_cn_fields;
                    cathode_cn_fields.reserve(state.cathode_particles.size());

                    for (auto &p : state.cathode_particles)
                    {
                        cathode_cn_fields.push_back(p.Cn_gf.get());
                    }

                    Utils::SaveSimulationSnapshotMulti(t, outdir, geometry, domain_parameters,
                        cathode_cn_fields, state.cathode_out, 5000);
                }

                t++;
            }

        
        
        // else
        // {
        //     RunFullCellSimulation(state, geometry, domain_parameters, bc, adjust, outdir, cfg);
        // }
        


        }
    }
    
    if (mfem::Mpi::WorldRank() == 0) { std::cout << "Simulation complete.\n"; }

    // End timing and output the total program execution time
    auto program_end = high_resolution_clock::now();
    if (mfem::Mpi::WorldRank() == 0) {std::cout << "Total Program Time: " 
            << duration_cast<seconds>(program_end - program_start).count() 
            << " seconds" << std::endl;}

*/
    }
    
    // Finalize HYPRE processing
    mfem::Hypre::Finalize();

    // Finalize MPI processing
    mfem::Mpi::Finalize();

    return 0;
}
