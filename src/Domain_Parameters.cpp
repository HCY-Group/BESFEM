#include "../include/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/readtiff.h"
#include "mfem.hpp"
#include <tiffio.h>
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

static inline void GlobalMinMax(const mfem::ParGridFunction& gf,
                                double& gmin, double& gmax,
                                MPI_Comm comm = MPI_COMM_WORLD)
{
    double lmin =  std::numeric_limits<double>::infinity();
    double lmax = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < gf.Size(); ++i) {
        const double v = gf(i);
        if (v < lmin) lmin = v;
        if (v > lmax) lmax = v;
    }
    MPI_Allreduce(&lmin, &gmin, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&lmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, comm);
}

double gTrgI = 0.0;

Domain_Parameters::Domain_Parameters(Initialize_Geometry &geo)
    : geometry(geo), nV(geo.nV), nE(geo.nE), nC(geo.nC), 
    dsF(geo.dsF   ? geo.dsF.get()   : nullptr),
    dsF_A(geo.dsF_A ? geo.dsF_A.get() : nullptr),
    dsF_C(geo.dsF_C ? geo.dsF_C.get() : nullptr),
    pmesh(geo.parallelMesh.get()), fespace(geo.parfespace)
{}

// Destructor
Domain_Parameters::~Domain_Parameters() {}

void Domain_Parameters::SetupDomainParameters(const char* mesh_type){

    InitializeGridFunctions();
    InterpolateDomainParameters(mesh_type);
    CalculatePhasePotentialsAndTargetCurrent();

    psi->SaveAsOne("psi");
    pse->SaveAsOne("pse");
    AvP->SaveAsOne("AvP");
    AvE->SaveAsOne("AvE");
    AvB->SaveAsOne("AvB");

    ps1->SaveAsOne("ps1");
    ps2->SaveAsOne("ps2");
    ps3->SaveAsOne("ps3");

    AvP_1->SaveAsOne("AvP_1");
    AvP_2->SaveAsOne("AvP_2");
    AvP_3->SaveAsOne("AvP_3");

    AvP_E_1->SaveAsOne("AvP_E_1");
    AvP_E_2->SaveAsOne("AvP_E_2");
    AvP_E_3->SaveAsOne("AvP_E_3");


    PrintInfo();
}

void Domain_Parameters::InitializeGridFunctions() {

    if (!fespace) {
        throw std::runtime_error("Finite element space is not initialized.");
    }

    psi = std::make_unique<mfem::ParGridFunction>(fespace.get());
    pse = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_0 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvB = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvE = std::make_unique<mfem::ParGridFunction>(fespace.get());

    AvP_1_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_2_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_3_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_1_2_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());

    AvP_E_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_E_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_E_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());

    AvP_all_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_all_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP_all_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());

    ps1   = std::make_unique<mfem::ParGridFunction>(fespace.get());
    ps2   = std::make_unique<mfem::ParGridFunction>(fespace.get());
    ps3   = std::make_unique<mfem::ParGridFunction>(fespace.get());

    Weight_all_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Weight_all_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Weight_all_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());

    Weight_1_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Weight_2_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Weight_3_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());

    Weight_E_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Weight_E_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Weight_E_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());

    psi_1_2 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    psi_2_3 = std::make_unique<mfem::ParGridFunction>(fespace.get());
    psi_3_1 = std::make_unique<mfem::ParGridFunction>(fespace.get());

    denom = std::make_unique<mfem::ParGridFunction>(fespace.get());

    const bool full = (dsF_A != nullptr) && (dsF_C != nullptr);
    if (full) {
        psA = std::make_unique<mfem::ParGridFunction>(fespace.get());
        psC = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvA = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvC = std::make_unique<mfem::ParGridFunction>(fespace.get());
    }

}

void Domain_Parameters::InterpolateDomainParameters(const char* mesh_type) {

    if (mesh_type == nullptr) {
        std::cerr << "Error: Mesh type not specified. Use -t option to specify mesh type (ml for matlab, v for voxel)." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    const bool full = (dsF_A != nullptr) && (dsF_C != nullptr);
    const bool is_root = (mfem::Mpi::WorldRank() == 0);
    const int dimension = pmesh->Dimension();

    if (!full) {
        // ------- HALF CELL (unchanged, but read from best available dsF) -------
        const mfem::ParGridFunction* g =
              dsF_A ? dsF_A
            : dsF_C ? dsF_C
            : dsF;

        if (!g) mfem::mfem_error("HALF mode: no active distance field available.");

        nV = pmesh->GetNV();

        if (strcmp(mesh_type, "ml") == 0) {
            for (int vi = 0; vi < nV; vi++) {
                (*psi)(vi) = 0.5 * (1.0 + tanh((*g)(vi) / (Constants::zeta * Constants::dh))); // matlab
                (*AvP)(vi) = -(pow(tanh((*g)(vi) / (Constants::zeta * Constants::dh)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // matlab

                (*pse)(vi) = 1.0 - (*psi)(vi);

                if ((*psi)(vi) < 0) { (*psi)(vi) = 0; }
                if ((*psi)(vi) > 1) { (*psi)(vi) = 1; }

                if ((*pse)(vi) < 0) { (*pse)(vi) = 0; }
                if ((*pse)(vi) > 1) { (*pse)(vi) = 1; }

                (*psi)(vi) += 1.0e-6; // Avoid zero values
                (*pse)(vi) += 1.0e-6; // Avoid zero values
            }
        }
        
        // if (strcmp(mesh_type, "v") == 0) {

        //     *psi = *geometry.MaskFilter;
        //     *pse = *geometry.MaskFilterPse;

        //     for (int i = 0; i < psi->Size(); i++)
        //     {

        //         if ((*psi)(i) < 0) { (*psi)(i) = 0; }
        //         if ((*psi)(i) > 1) { (*psi)(i) = 1; }

        //         if ((*pse)(i) < 0) { (*pse)(i) = 0; }
        //         if ((*pse)(i) > 1) { (*pse)(i) = 1; }

        //         (*psi)(i) += 1.0e-6; // Avoid zero values
        //         (*pse)(i) += 1.0e-6; // Avoid zero values

        //     }
        // }

        if (strcmp(mesh_type, "v") == 0) {

            MFEM_VERIFY(geometry.Vox, "geometry.Vox is not initialized.");

            *ps1 = *geometry.MaskFilter1;
            *ps2 = *geometry.MaskFilter2;
            *ps3 = *geometry.MaskFilter3;
            *pse = *geometry.MaskFilterPse;

            *psi = 0.0;
            *psi += *ps1;
            *psi += *ps2;
            *psi += *ps3;

            for (int i = 0; i < psi->Size(); i++)
            {
                if ((*psi)(i) < 0) { (*psi)(i) = 0.0; }
                if ((*psi)(i) > 1) { (*psi)(i) = 1.0; }

                if ((*pse)(i) < 0) { (*pse)(i) = 0.0; }
                if ((*pse)(i) > 1) { (*pse)(i) = 1.0; }

                if ((*ps1)(i) < 0) { (*ps1)(i) = 0.0; }
                if ((*ps1)(i) > 1) { (*ps1)(i) = 1.0; }

                if ((*ps2)(i) < 0) { (*ps2)(i) = 0.0; }
                if ((*ps2)(i) > 1) { (*ps2)(i) = 1.0; }

                if ((*ps3)(i) < 0) { (*ps3)(i) = 0.0; }
                if ((*ps3)(i) > 1) { (*ps3)(i) = 1.0; }

                (*psi)(i) += 1.0e-6;
                (*pse)(i) += 1.0e-6;
                (*ps1)(i) += 1.0e-6;
                (*ps2)(i) += 1.0e-6;
                (*ps3)(i) += 1.0e-6;
            }
        }

        // -------------------------------------------------
        // MULTI PARTICLE AvP
        // -------------------------------------------------

        if (strcmp(mesh_type, "v") == 0) 
        {

            auto ComputeGradMagnitude = [&](const mfem::ParGridFunction &phase_in,
                                            mfem::ParGridFunction &AvP_out)
            {
                const int dim = pmesh->Dimension();
                mfem::ParGridFunction dphase(fespace.get());

                AvP_out = 0.0;

                for (int d = 0; d < dim; ++d)
                {
                    dphase = 0.0;

                    // GetDerivative is non-const in practice, so copy if needed
                    mfem::ParGridFunction phase_tmp(phase_in);
                    phase_tmp.GetDerivative(1, d, dphase);

                    for (int vi = 0; vi < nV; ++vi)
                    {
                        const double v = dphase(vi);
                        AvP_out(vi) += v * v;
                    }
                }

                for (int vi = 0; vi < nV; ++vi)
                {
                    AvP_out(vi) = std::sqrt(AvP_out(vi));
                }
            };

            ComputeGradMagnitude(*psi, *AvP);      // total solid
            ComputeGradMagnitude(*ps1, *AvP_1);
            ComputeGradMagnitude(*ps2, *AvP_2);
            ComputeGradMagnitude(*ps3, *AvP_3);


            auto BuildPairInterface = [&](mfem::ParGridFunction &out,
                                const mfem::ParGridFunction &psa,
                                const mfem::ParGridFunction &psb,
                                const mfem::ParGridFunction &AvPa,
                                const mfem::ParGridFunction &AvPb)
            {
                out = psa;      // psa
                out *= AvPb;    // psa * AvPb

                mfem::ParGridFunction tmp(fespace.get());
                tmp = psb;      // psb
                tmp *= AvPa;    // psb * AvPa
                out += tmp;     // psa*AvPb + psb*AvPa

                mfem::ParGridFunction overlap(fespace.get());
                overlap = psa;
                overlap *= psb; // psa * psb

                out *= overlap;
                out *= 4.0;

                for (int vi = 0; vi < out.Size(); ++vi)
                {
                    if (out(vi) > 9000.0)
                    {
                        out(vi) = 1.4e4;
                    }
                }
            };

            BuildPairInterface(*AvP_1_2, *ps1, *ps2, *AvP_1, *AvP_2);
            BuildPairInterface(*AvP_2_3, *ps2, *ps3, *AvP_2, *AvP_3);
            BuildPairInterface(*AvP_3_1, *ps3, *ps1, *AvP_3, *AvP_1);

            auto BuildElectrolyteInterface = [&](mfem::ParGridFunction &out,
                                     const mfem::ParGridFunction &pse_in,
                                     const mfem::ParGridFunction &AvPk)
            {
                out = pse_in;
                out *= AvPk;
            };

            BuildElectrolyteInterface(*AvP_E_1, *pse, *AvP_1);
            BuildElectrolyteInterface(*AvP_E_2, *pse, *AvP_2);
            BuildElectrolyteInterface(*AvP_E_3, *pse, *AvP_3);

            *AvP_1_2_3 = 0.0;
            *AvP_1_2_3 += *AvP_1_2;
            *AvP_1_2_3 += *AvP_2_3;
            *AvP_1_2_3 += *AvP_3_1;

            *AvP_all_1 = 0.0;
            *AvP_all_1 += *AvP_1_2;
            *AvP_all_1 += *AvP_3_1;

            *AvP_all_2 = 0.0;
            *AvP_all_2 += *AvP_1_2;
            *AvP_all_2 += *AvP_2_3;

            *AvP_all_3 = 0.0;
            *AvP_all_3 += *AvP_3_1;
            *AvP_all_3 += *AvP_2_3;

            *denom = 0.0;
            *denom += *AvP_1_2_3;
            *denom += *AvP_E_1;
            *denom += *AvP_E_2;
            *denom += *AvP_E_3;

            *psi_1_2 = *ps1;
            *psi_1_2 += *ps2;
            for (int vi = 0; vi < nV; ++vi)
            {
                if ((*psi_1_2)(vi) > 1.0) { (*psi_1_2)(vi) = 1.0; }
            }

            *psi_2_3 = *ps2;
            *psi_2_3 += *ps3;
            for (int vi = 0; vi < nV; ++vi)
            {
                if ((*psi_2_3)(vi) > 1.0) { (*psi_2_3)(vi) = 1.0; }
            }

            *psi_3_1 = *ps3;
            *psi_3_1 += *ps1;
            for (int vi = 0; vi < nV; ++vi)
            {
                if ((*psi_3_1)(vi) > 1.0) { (*psi_3_1)(vi) = 1.0; }
            }

            auto ComputeWeight = [&](mfem::ParGridFunction &weight_out,
                         const mfem::ParGridFunction &num_in,
                         const mfem::ParGridFunction *mask_in = nullptr)
            {
                weight_out = 0.0;

                const double beta = 0.8;
                const double eps  = 1e-30;

                for (int vi = 0; vi < nV; ++vi)
                {
                    const double num = num_in(vi);
                    const double den = (*denom)(vi);

                    double ratio = 0.0;
                    if (den > eps)
                    {
                        ratio = num / den;
                        if (ratio < 0.0) ratio = 0.0;
                    }

                    weight_out(vi) = std::pow(ratio, beta);
                }

                if (mask_in) { weight_out *= *mask_in; }
            };

            ComputeWeight(*Weight_E_1, *AvP_E_1);
            ComputeWeight(*Weight_E_2, *AvP_E_2);
            ComputeWeight(*Weight_E_3, *AvP_E_3);

            ComputeWeight(*Weight_all_1, *AvP_all_1, psi.get());
            ComputeWeight(*Weight_all_2, *AvP_all_2, psi.get());
            ComputeWeight(*Weight_all_3, *AvP_all_3, psi.get());

            ComputeWeight(*Weight_1_2, *AvP_1_2, psi_1_2.get());
            ComputeWeight(*Weight_2_3, *AvP_2_3, psi_2_3.get());
            ComputeWeight(*Weight_3_1, *AvP_3_1, psi_3_1.get());

        }



        // // =====================================================
        // //  Calculating AvP for TIF Voxel
        // // =====================================================

        // if (strcmp(mesh_type, "v") == 0)
        // {
        //     const int dim = pmesh->Dimension(); 
        //     mfem::ParGridFunction dpsi(fespace.get());

        //     (*AvP) = 0.0;
        //     for (int d = 0; d < dim; d++)
        //     {
        //         dpsi = 0.0;
        //         psi->GetDerivative(1, d, dpsi); 

        //         // compound squares: AvP += (dpsi)^2
        //         for (int vi = 0; vi < nV; vi++)
        //         {
        //             const double v = dpsi(vi);
        //             (*AvP)(vi) += v * v;
        //         }
        //     }

        //     // sqrt to get magnitude
        //     for (int vi = 0; vi < nV; vi++)
        //     {
        //         (*AvP)(vi) = std::sqrt((*AvP)(vi));
        //     }

            
        //     // AvE
        //     mfem::ParGridFunction dpse(fespace.get());

        //     (*AvE) = 0.0;
        //     for (int d = 0; d < dim; d++)
        //     {
        //         dpse = 0.0;
        //         pse->GetDerivative(1, d, dpse); 

        //         // compound squares: AvE += (dpse)^2
        //         for (int vi = 0; vi < nV; vi++)
        //         {
        //             const double v = dpse(vi);
        //             (*AvE)(vi) += v * v;
        //         }
        //     }

        //     // sqrt to get magnitude
        //     for (int vi = 0; vi < nV; vi++)
        //     {
        //         (*AvE)(vi) = std::sqrt((*AvE)(vi));
        //     }
            
        // }

        // // =====================================================
        // //  End of Calculating AvP
        // // =====================================================

        // ---- GLOBAL checks for psi -------------------------------------------
        double psi_min = 0.0, psi_max = 0.0;
        GlobalMinMax(*psi, psi_min, psi_max);

        // Basic bounds check
        if (mfem::Mpi::WorldRank() == 0) {std::cout << "[Psi Check] min = " << psi_min 
                << ", max = " << psi_max << " (expected min = 1e-06, max = 1)" << std::endl;}
    

        if (psi_min < 0.0 || psi_max > 1.0 + 1e-6) {
            std::cerr << "[Psi Check] ERROR: psi values out of [0,1]!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (psi_min > 0.1) {
            std::cerr << "[Psi Check] ERROR: psi_min not near 0." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (psi_max < 0.9) {
            std::cerr << "[Psi Check] ERROR: psi_max not near 1." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // =====================================================
        //  Calculating AvB
        // =====================================================
        
        AvB = std::make_unique<mfem::ParGridFunction>(*AvP);

        if (strcmp(mesh_type, "v") == 0) {

            *AvB *= *AvE;

            for (int vi = 0; vi < nV; vi++)
            {
                (*AvB)(vi) = std::sqrt((*AvB)(vi));
            }

            for (int vi = 0; vi < psi->Size(); vi++)
            {
                const double psi_v = (*psi)(vi);
                const double pse_v = (*pse)(vi);

                const bool in_overlap =
                    (psi_v > 0.05 && psi_v < 0.95) &&
                    (pse_v > 0.05 && pse_v < 0.95);

                if (!in_overlap) { (*AvB)(vi) = 0.0; }
            }

        }

    } else {
        // ------- FULL CELL -------

        for (int vi = 0; vi < nV; vi++) {
            if (strcmp(mesh_type, "ml") == 0) {
                (*psA)(vi) = 0.5 * (1.0 + tanh((*dsF_A)(vi) / (Constants::zeta * Constants::dh))); // matlab
                (*AvA)(vi) = -(pow(tanh((*dsF_A)(vi) / (Constants::zeta * Constants::dh)), 2) - 1.0) / (Constants::zeta * Constants::dh); // matlab
                (*psC)(vi) = 0.5 * (1.0 + tanh((*dsF_C)(vi) / (Constants::zeta * Constants::dh))); // matlab
                (*AvC)(vi) = -(pow(tanh((*dsF_C)(vi) / (Constants::zeta * Constants::dh)), 2) - 1.0) / (Constants::zeta * Constants::dh); // matlab

            } else if (strcmp(mesh_type, "v") == 0) {
                (*psA)(vi) = 0.5 * (1.0 + tanh((*dsF_A)(vi))); // voxel
                (*AvA)(vi) = -(pow(tanh((*dsF_A)(vi)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // voxel
                (*psC)(vi) = 0.5 * (1.0 + tanh((*dsF_C)(vi))); // voxel
                (*AvC)(vi) = -(pow(tanh((*dsF_C)(vi)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // voxel

            }

            (*pse)(vi) = 1.0 - (*psA)(vi) - (*psC)(vi);

            if ((*psA)(vi) < Constants::eps) { (*psA)(vi) = Constants::eps; }
            if ((*pse)(vi) < Constants::eps) { (*pse)(vi) = Constants::eps; }
            if ((*psC)(vi) < Constants::eps) { (*psC)(vi) = Constants::eps; }

        }

            *psi = 0.0;
            *psi += *psA;
            *psi += *psC;

            double psA_min = 0, psA_max = 0, psC_min = 0, psC_max = 0;
            GlobalMinMax(*psA, psA_min, psA_max);
            GlobalMinMax(*psC, psC_min, psC_max);

            // Basic bounds check
            std::cout << "[PsA Check] min = " << psA_min 
                    << ", max = " << psA_max << " (expected min = 1e-06, max = 1)" << std::endl;

            if (psA_min < 0.0 || psA_max > 1.0 + 1e-6) {
                std::cerr << "[PsA Check] ERROR: psA values out of [0,1]!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (psA_min > 2e-6) {
                std::cerr << "[PsA Check] ERROR: psA_min not near 0." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (psA_max < 0.9) {
                std::cerr << "[PsA Check] ERROR: psA_max not near 1." << std::endl;
                std::exit(EXIT_FAILURE);
            }

            std::cout << "[psC Check] min = " << psC_min 
                << ", max = " << psC_max << " (expected min = 1e-06, max = 1)" << std::endl;

            if (psC_min < 0.0 || psC_max > 1.0 + 1e-6) {
                std::cerr << "[psC Check] ERROR: psC values out of [0,1]!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (psC_min > 2e-6) {
                std::cerr << "[psC Check] ERROR: psC_min not near 0." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (psC_max < 0.9) {
                std::cerr << "[psC Check] ERROR: psC_max not near 1." << std::endl;
                std::exit(EXIT_FAILURE);
            }

        
        AvB = std::make_unique<mfem::ParGridFunction>(*AvA);
        
        for (int vi = 0; vi < nV; vi++) {
            if ((*AvA)(vi) * Constants::dh < Constants::thres) { (*AvA)(vi) = 0.0; }
            if ((*AvC)(vi) * Constants::dh < Constants::thres) { (*AvC)(vi) = 0.0; }
            if ((*AvB)(vi) * Constants::dh < 1.0e-5) { (*AvB)(vi) = 0.0; }
        }

    }
}

void Domain_Parameters::CalculateTotals(const mfem::ParGridFunction& grid_function, const mfem::Vector& element_volumes, double& local_total, double& global_total) {
    local_total = 0.0;

    // Average value for each element
    mfem::Vector element_avg_values(nE);

    for (int ei = 0; ei < pmesh->GetNE(); ei++) {
        mfem::Array<double> vertex_values(nC);
        grid_function.GetNodalValues(ei, vertex_values);

        // Compute the average value over all corners of the element
        double avg_value = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            avg_value += vertex_values[vt];
        }
        avg_value /= nC;

        element_avg_values(ei) = avg_value;

        // Accumulate the weighted value for the current element
        local_total += avg_value * element_volumes(ei);
    }

    // Perform global MPI reduction to sum values across all processes
    MPI_Allreduce(&local_total, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void Domain_Parameters::CalculateTotalPhaseField(const mfem::ParGridFunction& grid_function, double& total, double& global_total) {

    EVol.SetSize(nE);
	for (int ei = 0; ei < nE; ei++){
        EVol(ei) = pmesh->GetElementVolume(ei);
        // EVol(ei) /= Constants::dh * Constants::dh; // in 2D
	}

    // Call the general `CalculateTotals` method for any field
    CalculateTotals(grid_function, EVol, total, global_total);

}

void Domain_Parameters::CalculatePhasePotentialsAndTargetCurrent() {

    const bool full = (dsF_A != nullptr) && (dsF_C != nullptr);

    // Half Cell : use psi
    if (!full) {
        // Calculate totals for Psi and Pse fields
        CalculateTotalPhaseField(*psi, tPsi, gtPsi);
        CalculateTotalPhaseField(*pse, tPse, gtPse);
        CalculateTargetCurrent(tPsi, gTrgI);

        CalculateTotalPhaseField(*ps1, tPsi1, gtPsi1);
        CalculateTotalPhaseField(*ps2, tPsi2, gtPsi2);
        CalculateTotalPhaseField(*ps3, tPsi3, gtPsi3);

        CalculateTargetCurrent(tPsi1, gTrg1);
        CalculateTargetCurrent(tPsi2, gTrg2);
        CalculateTargetCurrent(tPsi3, gTrg3);
    }

    // Full Cell : use psA, psC
    else {
        // Calculate totals for PsA, PsC, and Pse fields
        CalculateTotalPhaseField(*psA, tPsA, gtPsA);
        CalculateTotalPhaseField(*psC, tPsC, gtPsC);
        CalculateTotalPhaseField(*pse, tPse, gtPse);
        CalculateTargetCurrent(tPsC, gTrgI);
    }

}

// void Domain_Parameters::CalculateTargetCurrent(double total_psi) {

//     // Compute target current based on total Psi, rho, Cr, and constants
//     trgI = total_psi * Constants::rho_C * (0.95 - 0.3) / (3600.0 / Constants::Cr); // bounds of cathode 

//     // Perform global MPI reduction to get the total target current
//     MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// }

void Domain_Parameters::CalculateTargetCurrent(double total_psi, double &global_total) {

    // Compute target current based on total Psi, rho, Cr, and constants
    trgI = total_psi * Constants::rho_C * (0.95 - 0.3) / (3600.0 / Constants::Cr); // bounds of cathode 

    // voxel
    // trgI = (total_psi * (pow(Constants::dh, 2)))  * Constants::rho_C * (0.95 - 0.3) / (3600.0 / Constants::Cr); // bounds of cathode 

    // Perform global MPI reduction to get the total target current
    MPI_Allreduce(&trgI, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void Domain_Parameters::PrintInfo() {

    const bool full = (dsF_A != nullptr) && (dsF_C != nullptr);

    if (mfem::Mpi::WorldRank() == 0) // only print on rank 0
    if (!full)
    {
        cout << "Total Psi: " << gtPsi << endl;
        cout << "Total Pse: " << gtPse << endl;
        cout << "Target Current: " << gTrgI << endl;
    }
    else
    {
        cout << "Total PsA: " << gtPsA << endl;
        cout << "Total PsC: " << gtPsC << endl;
        cout << "Total Pse: " << gtPse << endl;
        cout << "Target Current: " << gTrgI << endl;
    }
}

