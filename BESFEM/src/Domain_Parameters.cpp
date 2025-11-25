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

static inline bool AnyRankTrue(bool local, MPI_Comm comm = MPI_COMM_WORLD)
{
    int l = local ? 1 : 0, g = 0;
    MPI_Allreduce(&l, &g, 1, MPI_INT, MPI_MAX, comm);
    return g != 0;
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

    PrintInfo();
}

void Domain_Parameters::InitializeGridFunctions() {

    if (!fespace) {
        throw std::runtime_error("Finite element space is not initialized.");
    }

    psi = std::make_unique<mfem::ParGridFunction>(fespace.get());
    pse = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvB = std::make_unique<mfem::ParGridFunction>(fespace.get());

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


    if (!full) {
        // ------- HALF CELL (unchanged, but read from best available dsF) -------
        const mfem::ParGridFunction* g =
              dsF_A ? dsF_A
            : dsF_C ? dsF_C
            : dsF;

        if (!g) mfem::mfem_error("HALF mode: no active distance field available.");


        if (strcmp(mesh_type, "v") == 0) {
            psi->ProjectGridFunction(*g);
            *psi -= 0.5;
        }

        for (int vi = 0; vi < nV; vi++) {
            if (strcmp(mesh_type, "ml") == 0) {
                (*psi)(vi) = 0.5 * (1.0 + tanh((*g)(vi) / (Constants::zeta * Constants::dh))); // matlab
                (*AvP)(vi) = -(pow(tanh((*g)(vi) / (Constants::zeta * Constants::dh)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // matlab

            } else if (strcmp(mesh_type, "v") == 0) {
                (*psi)(vi) = 0.5 * (1.0 + tanh((*g)(vi))); // voxel
                (*AvP)(vi) = -(pow(tanh((*g)(vi)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // voxel

            } 

            (*pse)(vi) = 1.0 - (*psi)(vi);

            if ((*psi)(vi) < 0) { (*psi)(vi) = 0; }
            if ((*psi)(vi) > 1) { (*psi)(vi) = 1; }

            if ((*pse)(vi) < 0) { (*pse)(vi) = 0; }
            if ((*pse)(vi) > 1) { (*pse)(vi) = 1; }

            (*psi)(vi) += 1.0e-6; // Avoid zero values
            (*pse)(vi) += 1.0e-6; // Avoid zero values

        }

            // ---- GLOBAL checks for psi -------------------------------------------
            double psi_min = 0.0, psi_max = 0.0;
            GlobalMinMax(*psi, psi_min, psi_max);

            // Basic bounds check
            std::cout << "[Psi Check] min = " << psi_min 
                    << ", max = " << psi_max << " (expected min = 1e-06, max = 1)" << std::endl;
        

            if (psi_min < 0.0 || psi_max > 1.0 + 1e-6) {
                std::cerr << "[Psi Check] ERROR: psi values out of [0,1]!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (psi_min > 2e-6) {
                std::cerr << "[Psi Check] ERROR: psi_min not near 0." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (psi_max < 1.0) {
                std::cerr << "[Psi Check] ERROR: psi_max not near 1." << std::endl;
                std::exit(EXIT_FAILURE);
            }

        
        AvB = std::make_unique<mfem::ParGridFunction>(*AvP);


        for (int vi = 0; vi < nV; vi++) {
            if ((*AvP)(vi) * Constants::dh < 1.0e-2) { (*AvP)(vi) = 0.0; }
            if ((*AvB)(vi) * Constants::dh < 1.0e-6) { (*AvB)(vi) = 0.0; }
        }
        

        mfem::GridFunctionCoefficient AvP_coeff(AvP.get());
        AvP->ProjectCoefficient(AvP_coeff);

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
        
        // mfem::GridFunctionCoefficient AvA_coeff(AvA.get());
        // AvA->ProjectCoefficient(AvA_coeff);

        // mfem::GridFunctionCoefficient AvC_coeff(AvC.get());
        // AvC->ProjectCoefficient(AvC_coeff);

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
        CalculateTargetCurrent(tPsi);
    }

    // Full Cell : use psA, psC
    else {
        // Calculate totals for PsA, PsC, and Pse fields
        CalculateTotalPhaseField(*psA, tPsA, gtPsA);
        CalculateTotalPhaseField(*psC, tPsC, gtPsC);
        CalculateTotalPhaseField(*pse, tPse, gtPse);
        CalculateTargetCurrent(tPsC);
    }

}

void Domain_Parameters::CalculateTargetCurrent(double total_psi) {

    // Compute target current based on total Psi, rho, Cr, and constants
    trgI = total_psi * Constants::rho_C * (0.95 - 0.3) / (3600.0 / Constants::Cr); // bounds of cathode 

    // Perform global MPI reduction to get the total target current
    MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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

