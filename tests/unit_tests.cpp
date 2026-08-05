// Unit tests and benchmark tests
// compile with make test
// run with: ./unit_tests

// Uses Catch2 (version 2.13)

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
// -probably at least one test for each source file...
// -do the domain parameters and interfaces make sense?
    // AvP_Pairs? Weights?
// -should do some 3d input files
// -should check if tests work in parallel with multiple cpus or gpus
// -should check if AMR works
// -maybe do some more simple multi-particle tests?

#include "mfem.hpp"
#include "mpi.h"
#include "../include/BESFEM_All.hpp"
#include "catch.hpp"

#include <chrono>
#include <iostream>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <vector>


TEST_CASE("Domain Parameters and Interfaces", "[psi]") {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto config_file = GENERATE("../tests/test_run_config_cathode.txt",
                                "../tests/test_run_config_anode.txt",
                                "../tests/test_run_config_cathode_multipart.txt");
    SimulationConfig cfg = SimConfigFromFileName(config_file);
    ValidateConfig(cfg, 0, nullptr);

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

    int np = domain_parameters.ps.size();

    mfem::ParGridFunction diff(*domain_parameters.psi);
   
    // sum of particle order parameters should be solid phase
    // sum(ps[k]) = psi
    diff = 0.0;
    for (int k=0; k<np; k++) {
        diff += *domain_parameters.ps[k];
        std::cout << "max diff: " << diff.Max() << std::endl;
    }
    diff -= *domain_parameters.psi;
    diff.SaveAsOne("psi_diff");
    CHECK( diff.Norml2() < 1.0e-3 );

    // range of values for psi and ps[k] and pse should be 1e-6 and 1.0
    CHECK( abs(domain_parameters.psi->Max() - 1.0) < 1.0e-3 );
    CHECK( abs(domain_parameters.pse->Max() - 1.0) < 1.0e-3 );
    for (int k=0; k<np; k++) {
        CHECK( abs(domain_parameters.ps[k]->Max() - 1.0) < 1.0e-3 );
    }
    CHECK( domain_parameters.psi->Min() > 0.0 );
    CHECK( domain_parameters.pse->Min() > 0.0 );
    for (int k=0; k<np; k++) {
        CHECK( domain_parameters.psi->Min() > 0.0 );
    }
    CHECK( domain_parameters.psi->Min() < 1.0e-3 );
    CHECK( domain_parameters.pse->Min() < 1.0e-3 );
    for (int k=0; k<np; k++) {
        CHECK( domain_parameters.psi->Min() < 1.0e-3 );
    }

    // magnitude of AvE and AvP and AvEs[k] and AvPs[k] should be approximately the same
    double diff_flt;
    diff_flt = domain_parameters.AvP->Max() - domain_parameters.AvE->Max();
    diff_flt /= domain_parameters.AvP->Max(); //normalize by max(grad(psi))
    CHECK( diff_flt < 0.05 );
    for (int k=0; k<np; k++) {
        diff_flt = domain_parameters.AvP->Max() - domain_parameters.AvPs[k]->Max();
        diff_flt /= domain_parameters.AvP->Max(); //normalize by max(grad(psi))
        CHECK( diff_flt < 0.05 );
        diff_flt = domain_parameters.AvP->Max() - domain_parameters.AvEs[k]->Max();
        diff_flt /= domain_parameters.AvP->Max(); //normalize by max(grad(psi))
        CHECK( diff_flt < 0.05 );
    }

/*
    for (int j=0; j<np; j++) {
        string weight_e_name = "Weight_Es" + std::to_string(j);
        domain_parameters.WeightEs[j]->SaveAsOne(weight_e_name.c_str());
        for (int k=j+1; k<np; k++) {
            string AvP_Pairs_name = "AvP_Pairs" + std::to_string(j) + std::to_string(k);
            domain_parameters.AvP_Pairs[j][k]->SaveAsOne(AvP_Pairs_name.c_str());
            string weight_pairs_name = "Weight_Pairs" + std::to_string(j) + std::to_string(k);
            domain_parameters.WeightPairs[j][k]->SaveAsOne(weight_pairs_name.c_str());
        }
    }
            




    // difference of particle interfaces and particle-particle interfaces should be bulk solid interface:
        // AvPs[k]-AvP_Pairs[i][k] = AvP?
    diff = 0.0;
    domain_parameters.AvP->SaveAsOne("AvP");
    for (int k=0; k<np; k++) {
        string AvEs_name = "AvEs" + std::to_string(k);
        string AvPs_name = "AvPs" + std::to_string(k);
        domain_parameters.AvEs[k]->SaveAsOne(AvEs_name.c_str());
        domain_parameters.AvPs[k]->SaveAsOne(AvPs_name.c_str());
        diff += *domain_parameters.AvPs[k];
        std::cout << "max diff: " << diff.Max() << std::endl;
    }
    diff.SaveAsOne("avps_sum");
    diff -= *domain_parameters.AvP;
    diff.SaveAsOne("avp_diff");
    CHECK( diff.Norml2() < 1.0e-3 );
 
    // sum of electrolyte-particle interfaces should be AvE: 
        // sum(AvEs[k]) = AvE
    diff = 0.0;
    domain_parameters.AvE->SaveAsOne("AvE");
    for (int k=0; k<np; k++) {
        string AvEs_name = "AvEs" + std::to_string(k);
        string AvPs_name = "AvPs" + std::to_string(k);
        domain_parameters.AvEs[k]->SaveAsOne(AvEs_name.c_str());
        domain_parameters.AvPs[k]->SaveAsOne(AvPs_name.c_str());
        diff += *domain_parameters.AvEs[k];
        std::cout << "max diff: " << diff.Max() << std::endl;
    }
    diff.SaveAsOne("aves_sum");
    diff -= *domain_parameters.AvE;
    diff.SaveAsOne("ave_diff");
    CHECK( diff.Norml2() < 1.0e-3 );
*/

}






TEST_CASE("UpdatePotential", "[phi]") {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    {

    auto config_file = GENERATE("../tests/test_run_config_cathode.txt",
                                "../tests/test_run_config_anode.txt");
    SimulationConfig cfg = SimConfigFromFileName(config_file);
    ValidateConfig(cfg, 0, nullptr);

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
        {
        std::vector<mfem::ParGridFunction*> cn_fields; // vector of pointers to cathode concentration fields
        std::vector<mfem::ParGridFunction*> psi_fields; // vector of pointers to cathode potential fields
        std::vector<sim::MaterialType> materials; // vector of cathode material types
                    
        int np;
        if (cfg.half_electrode == sim::Electrode::ANODE){
            np = static_cast<int>(state.anode_particles.size());
        } else {
            np = static_cast<int>(state.cathode_particles.size());
        }
 
        cn_fields.reserve(np); // pre-allocate memory
        psi_fields.reserve(np); // pre-allocate memory
        materials.reserve(np); // pre-allocate memory

        for (int j = 0; j < np; ++j)
        {
            psi_fields.push_back(domain_parameters.ps[j].get());
            if (cfg.half_electrode == sim::Electrode::ANODE){
                cn_fields.push_back(state.anode_particles[j].Cn_gf.get()); 
                materials.push_back(state.anode_particles[j].material);
            } else {
                cn_fields.push_back(state.cathode_particles[j].Cn_gf.get()); 
                materials.push_back(state.cathode_particles[j].material);
            }
        }
        if (cfg.half_electrode == sim::Electrode::ANODE){
            state.anode_potential->AssembleSystem(cn_fields, psi_fields, materials, *state.phA_gf);
        } else {
            state.cathode_potential->AssembleSystem(cn_fields, psi_fields, materials, *state.phC_gf);
        }
        state.electrolyte_potential->AssembleSystem(*state.CnE_gf, *domain_parameters.pse, *state.phE_gf);

        double globalerror_P = 1.0; // Error for particle potential
        double globalerror_E = 1.0; // Error for electrolyte potential

        mfem::ParGridFunction Rxn_p(*domain_parameters.AvP);
        mfem::ParGridFunction Rxn_e(*domain_parameters.AvE);

        if (cfg.half_electrode == sim::Electrode::ANODE){
            state.anode_potential->UpdatePotential(Rxn_p, *state.phA_gf, *domain_parameters.psi, globalerror_P);
        } else {
            state.cathode_potential->UpdatePotential(Rxn_p, *state.phC_gf, *domain_parameters.psi, globalerror_P);
        }
        state.electrolyte_potential->UpdatePotential(Rxn_e, *state.phE_gf, *domain_parameters.pse, globalerror_E);
        
        domain_parameters.psi->SaveAsOne("psi_test");
        domain_parameters.pse->SaveAsOne("pse_test");
        if (cfg.half_electrode == sim::Electrode::ANODE){
            state.phA_gf->SaveAsOne("phiA");
        } else {
            state.phC_gf->SaveAsOne("phiC");
        }
        state.phE_gf->SaveAsOne("phiE");
        //std::cout << "Faraday: " << Constants::Frd << std::endl;

        // ANALYTICAL SOLUTION
        Rxn_p.FESpace()->GetMesh()->EnsureNodes();
        mfem::GridFunction *nodes = Rxn_p.FESpace()->GetMesh()->GetNodes();

        mfem::ParGridFunction phP_an(Rxn_p);
        mfem::ParGridFunction phE_an(Rxn_p);
        mfem::ParGridFunction x(Rxn_p);
        mfem::ParGridFunction y(Rxn_p);
        for (int i=0; i<phP_an.Size(); i++)
        {
            x(i) = (*nodes)(2*i);
            y(i) = (*nodes)(2*i+1);
        }
        
        if (cfg.half_electrode == sim::Electrode::ANODE){
            std::cout << "kappa: " << state.anode_potential->GetConductivity().Min() << " " << 
                                      state.anode_potential->GetConductivity().Max() << std::endl;
            std::cout << "F/kap: " << Constants::Frd/state.anode_potential->GetConductivity().Max() << std::endl;
        } else {
            std::cout << "kappa: " << state.cathode_potential->GetConductivity().Min() << " " << 
                                      state.cathode_potential->GetConductivity().Max() << std::endl;
            std::cout << "F/kap: " << Constants::Frd/state.cathode_potential->GetConductivity().Max() << std::endl;
        }


        double slope_e = 1.0/state.electrolyte_potential->GetConductivity().Max(); 
        double intercept_e = state.electrolyte_potential->GetBoundaryVoltage();
        double slope_p;
        double intercept_p;
        if (cfg.half_electrode == sim::Electrode::ANODE){
            slope_p = Constants::Frd/state.anode_potential->GetConductivity().Max(); 
            intercept_p = state.anode_potential->GetBoundaryVoltage();
        } else {
            slope_p = Constants::Frd/state.cathode_potential->GetConductivity().Max(); 
            intercept_p = state.cathode_potential->GetBoundaryVoltage();
        }
        std::cout << "particle slope and intercept: " << slope_p << " " << intercept_p << std::endl;
        for (int i=0; i < phP_an.Size(); i++)
        {
          if (cfg.half_electrode == sim::Electrode::ANODE){
            phP_an(i) = slope_p*y(i) + intercept_p;
            phE_an(i) = -slope_e*(y.Max()-y(i)) + intercept_e;
          } else {
            phP_an(i) = slope_p*(y.Max()-y(i)) + intercept_p;
            phE_an(i) = -slope_e*y(i) + intercept_e;
          }
        }
        phP_an.SaveAsOne("phiP_an");

        //Multiply by psi to compare
        phP_an *= *domain_parameters.psi;
        if (cfg.half_electrode == sim::Electrode::ANODE){
            *state.phA_gf *= *domain_parameters.psi;
            state.phA_gf->SaveAsOne("phiA");
        } else {
            *state.phC_gf *= *domain_parameters.psi;
            state.phC_gf->SaveAsOne("phiC");
        }
        phP_an.SaveAsOne("phiP_an");

        phE_an *= *domain_parameters.pse;
        *state.phE_gf *= *domain_parameters.pse;
        state.phE_gf->SaveAsOne("phiE");
        phE_an.SaveAsOne("phiE_an");

        // compute difference, calculate l2 norm
        mfem::ParGridFunction diff_c(phE_an);
        if (cfg.half_electrode == sim::Electrode::ANODE){
            diff_c = *state.phA_gf;
        } else {
            diff_c = *state.phC_gf;
        }
        diff_c -= phP_an;
        diff_c /= phP_an;  //weight by magnitude of analytical solution
        diff_c *= *domain_parameters.psi; // only get errors in domain of interest
        std::cout << "L2 error cathode:     " << diff_c.Norml2() << std::endl;

        mfem::ParGridFunction diff_e(*state.phE_gf);
        diff_e = *state.phE_gf;
        diff_e -= phE_an;
        diff_e /= phE_an;  //weight by magnitude of analytical solution
        diff_e *= *domain_parameters.pse; // only get errors in domain of interest
        std::cout << "L2 error electrolyte: " << diff_e.Norml2() << std::endl;

        REQUIRE( diff_c.Norml2() < 0.05 );
        REQUIRE( diff_e.Norml2() < 0.05 );


        }

    }
}



TEST_CASE("UpdateConcentration", "[Cn]") {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    {

    auto config_file = GENERATE("../tests/test_run_config_anode.txt",
                                "../tests/test_run_config_cathode_NMC.txt");
    SimulationConfig cfg = SimConfigFromFileName(config_file);
    ValidateConfig(cfg, 0, nullptr);

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


        // =================================
        // UPDATE CONCENTRATION
        // =================================
        for (int t=0; t<10; t++) {
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
                Pairs(state, geometry, domain_parameters, j, pair_terms, np, 0);

                //state.anode_particles[j].concentration->UpdateConcentration(*state.anode_particles[j].Rx_src, *state.anode_particles[j].Cn_gf,
                //   *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
                state.anode_particles[j].concentration->UpdateConcentration(*domain_parameters.AvEs[j], *state.anode_particles[j].Cn_gf,
                    *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
        
                string CnP_name = "CnP" + std::to_string(j);
                state.anode_particles[j].Cn_gf->SaveAsOne(CnP_name.c_str());

            }

            //state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf,
            //    *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});
            state.electrolyte_concentration->UpdateConcentration(*domain_parameters.AvP, *state.CnE_gf,
                *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});
            state.CnE_gf->SaveAsOne("CnE");


        }
        else 
        {
            const int np = static_cast<int>(state.cathode_particles.size());
            std::vector<double> global_currents(np, 0.0);

            UpdateCathodePairChemicalPotentials(state, geometry, domain_parameters);



            *state.Rxn_gf = 0.0;
            double Rxn_const = 1e-7;
            for (int j = 0; j < np; ++j)
            {
                // constant reaction (no Butler Volmer)
                *state.cathode_particles[j].Rxn_gf = *domain_parameters.AvEs[j];
                *state.cathode_particles[j].Rxn_gf *= Rxn_const;

                *state.cathode_particles[j].Rx_src = *state.cathode_particles[j].Rxn_gf;
                *state.Rxn_gf += *state.cathode_particles[j].Rxn_gf;

                std::vector<ConcentrationBase::PairCoupling> pair_terms;
                Pairs(state, geometry, domain_parameters, j, pair_terms, np, 0);

                state.cathode_particles[j].concentration->UpdateConcentration(*state.cathode_particles[j].Rx_src, *state.cathode_particles[j].Cn_gf,
                   *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
                //state.cathode_particles[j].concentration->UpdateConcentration(*domain_parameters.AvEs[j], *state.cathode_particles[j].Cn_gf,
                //    *domain_parameters.ps[j], domain_parameters.gtPs[j], *domain_parameters.WeightEs[j], pair_terms);
        
                string CnP_name = "CnP" + std::to_string(j);
                state.cathode_particles[j].Cn_gf->SaveAsOne(CnP_name.c_str());

            }
            state.Rxn_gf->SaveAsOne("Rxn_test");

            state.electrolyte_concentration->UpdateConcentration(*state.Rxn_gf, *state.CnE_gf,
                *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});
            //state.electrolyte_concentration->UpdateConcentration(*domain_parameters.AvP, *state.CnE_gf,
            //    *domain_parameters.pse, domain_parameters.gtPse, *domain_parameters.pse, {});
            state.CnE_gf->SaveAsOne("CnE");


            // ====================================
            // ANALYTICAL SOLUTION
            // semi-infinite solid, with constant heat flux at boundary
            // a d^2T/dx^2 = dT/dt
            //
            // Wiśniewski, T.S. (2014). Transient Heat Conduction in Semi-infinite Solid with Specified Surface Heat Flux. 
            // In: Hetnarski, R.B. (eds) Encyclopedia of Thermal Stresses. Springer, Dordrecht. 
            // https://doi.org/10.1007/978-94-007-2739-7_413
            //
            // Charles H. Forsberg, Chapter 4 - Unsteady conduction,
            // Heat Transfer Principles and Applications, Academic Press, 2021, Pages 121-162,
            // https://doi.org/10.1016/B978-0-12-802296-2.00004-4.
            // 
            // T(x,t) = T0 + q/k[ sqrt(4at/pi)exp(-x^2/4at) - x erfc(x/2sqrt(at)) ]
            // q = heat flux
            // k = thermal conductivity
            // q/k = dT/dx at boundary (fourier's law) - this is our imposed boundary condition
            // a = thermal diffusivity = k/(rho*c)
            // T0 = initial temp
            // t = time
            // x = depth
            // ====================================

        state.CnE_gf->FESpace()->GetMesh()->EnsureNodes();
        mfem::GridFunction *nodes = state.CnE_gf->FESpace()->GetMesh()->GetNodes();

        mfem::ParGridFunction x(*state.CnE_gf);
        mfem::ParGridFunction y(*state.CnE_gf);
        for (int i=0; i<state.CnE_gf->Size(); i++)
        {
            x(i) = (*nodes)(2*i);
            y(i) = (*nodes)(2*i+1);
        }
            
            mfem::ParGridFunction CnE_an(*state.CnE_gf);
            double diff_e = state.electrolyte_concentration->GetDiffusivity().Max(); 

            CnE_an = cfg.init_CnE;
            double time_elapsed = (t+1)*cfg.dt;
            // boundary condition (=q/k)
            double B_n = -Rxn_const*Constants::t_minus;                
            const double pi = std::acos(-1.0);
            for (int i=0; i<CnE_an.Size(); i++) {
                double a = std::sqrt( 4*diff_e*time_elapsed/pi );
                double b = std::exp( -x(i)*x(i)/4/diff_e/time_elapsed );
                double c = x(i)*( 1.0-std::erf( x(i)/2.0/std::sqrt(diff_e*time_elapsed)  )  );

                
                CnE_an(i) += B_n * (a*b - c);
                //std::cout << "a: " << a << ", b: " << b << ", c: " << c << std::endl;
                std::cout << CnE_an(i) << std::endl;
            }
            std::cout << "diffusivity: " << diff_e << std::endl;
            std::cout << "Boundary: " << B_n << std::endl;
            CnE_an.SaveAsOne("CnE_an");


        }
        }

    }
}

















