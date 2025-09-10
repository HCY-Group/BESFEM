#include "../inputs/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/readtiff.h"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>


double gTrgI = 0.0;

Domain_Parameters::Domain_Parameters(Initialize_Geometry &geo)
    : geometry(geo), nV(geo.nV), nE(geo.nE), nC(geo.nC), dsF(geo.dsF.get()), pmesh(geo.parallelMesh.get()),
    fespace(geo.parfespace)
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

    psi = make_unique<mfem::ParGridFunction>(fespace.get());
    pse = make_unique<mfem::ParGridFunction>(fespace.get());
    AvP = make_unique<mfem::ParGridFunction>(fespace.get());
    AvB = make_unique<mfem::ParGridFunction>(fespace.get());

}

void Domain_Parameters::InterpolateDomainParameters(const char* mesh_type) {

    if (mesh_type == nullptr) {
        std::cerr << "Error: Mesh type not specified. Use -t option to specify mesh type (r for rectangle, c for circle, v for voxel)." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (strcmp(mesh_type, "v") == 0) {
        psi->ProjectGridFunction(*dsF);
        *psi -= 0.5;
    }

    for (int vi = 0; vi < nV; vi++) {
        if (strcmp(mesh_type, "r") == 0) {
            (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (Constants::zeta * Constants::dh))); // rectangle
            (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (Constants::zeta * Constants::dh)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // rectangle

        } else if (strcmp(mesh_type, "c") == 0) {
            (*psi)(vi) = 0.5 * (1.0 + tanh((11.0 * Constants::dh - (*dsF)(vi)) / Constants::zeta)); // circle
            (*AvP)(vi) = -(pow(tanh((11.0 * Constants::dh - (*dsF)(vi)) / Constants::zeta), 2) - 1.0) / (2 * Constants::zeta); // circle

        } else if (strcmp(mesh_type, "d") == 0) {
            (*psi)(vi) = 0.5 * (1.0 + tanh(((*dsF)(vi)) / (Constants::zeta))); // disk
            (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (Constants::zeta)), 2) - 1.0) / (2 * Constants::zeta); // disk

        } else if (strcmp(mesh_type, "v") == 0) {
            (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi))); // voxel
            (*AvP)(vi) = -(pow(tanh((*dsF)(vi)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // voxel

        } else {
            (*psi)(vi) = 0.5 * (1.0 + tanh((*psi)(vi))); // voxel
            (*AvP)(vi) = -(pow(tanh((*dsF)(vi)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh); // voxel
        }

        (*pse)(vi) = 1.0 - (*psi)(vi);

        if ((*psi)(vi) < 0) { (*psi)(vi) = 0; }
        if ((*psi)(vi) > 1) { (*psi)(vi) = 1; }

        if ((*pse)(vi) < 0) { (*pse)(vi) = 0; }
        if ((*pse)(vi) > 1) { (*pse)(vi) = 1; }

        (*psi)(vi) += 1.0e-6; // Avoid zero values
        (*pse)(vi) += 1.0e-6; // Avoid zero values

    }

        double psi_min = psi->Min();
        double psi_max = psi->Max();

        // Basic bounds check
        std::cout << "[Psi Check] min = " << psi_min 
                << ", max = " << psi_max << " (expected min = 1e-06, max = 1)" << std::endl;

        if (psi_min < 0.0 || psi_max > 1.0 + 1e-6) {
            std::cerr << "[Psi Check] ERROR: psi values out of [0,1]!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (psi_min > 2e-6) {
            std::cerr << "[Psi Check] WARNING: psi_min not near 0." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (psi_max < 1.0) {
            std::cerr << "[Psi Check] WARNING: psi_max not near 1." << std::endl;
            std::exit(EXIT_FAILURE);
        }

    
    AvB = std::make_unique<mfem::ParGridFunction>(*AvP);

    if(strcmp(mesh_type, "d") == 0) 
        for (int vi = 0; vi < nV; vi++) {
            if ((*AvP)(vi) < 1.0e-2) { (*AvP)(vi) = 0.0; }
            if ((*AvB)(vi) < 1.0e-6) { (*AvB)(vi) = 0.0; }
        }    
    else {
        for (int vi = 0; vi < nV; vi++) {
            if ((*AvP)(vi) * Constants::dh < 1.0e-3) { (*AvP)(vi) = 0.0; }
            if ((*AvB)(vi) * Constants::dh < 1.0e-6) { (*AvB)(vi) = 0.0; }
        }
    }
    


    if(strcmp(mesh_type, "d") == 0) {
        *AvP /= Constants::dh;
    }

    mfem::GridFunctionCoefficient AvP_coeff(AvP.get());
    AvP->ProjectCoefficient(AvP_coeff);

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

    // Calculate totals for Psi and Pse fields
    CalculateTotalPhaseField(*psi, tPsi, gtPsi);
    CalculateTotalPhaseField(*pse, tPse, gtPse);

    // Compute the target current using Psi
    CalculateTargetCurrent(tPsi);
}

void Domain_Parameters::CalculateTargetCurrent(double total_psi) {

    // Compute target current based on total Psi, rho, Cr, and constants
    trgI = total_psi * Constants::rho * (0.9 - 0.3) / (3600.0 / Constants::Cr);
    // trgI = total_psi * Constants::rho * (1.0 - 0.0) / (3600.0 / Constants::Cr);

    // Perform global MPI reduction to get the total target current
    MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void Domain_Parameters::PrintInfo() {
    if (mfem::Mpi::WorldRank() == 0) // only print on rank 0
    {
        cout << "Total Psi: " << gtPsi << endl;
        cout << "Total Pse: " << gtPse << endl;
        cout << "Target Current: " << gTrgI << endl;
    }
}
