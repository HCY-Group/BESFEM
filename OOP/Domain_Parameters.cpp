#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "Constants.hpp"
#include "readtiff.h"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

double gTrgI = 0.0;

// Constructor
Domain_Parameters::Domain_Parameters(Initialize_Geometry &geo) 
    : geometry(geo), nV(geo.nV), nE(geo.nE), nC(geo.nC), dsF(geo.dsF.get()), pmesh(geo.parallelMesh.get()),
    fespace(geo.parfespace)
{}

// Destructor
Domain_Parameters::~Domain_Parameters() {}

void Domain_Parameters::SetupDomainParameters(){
    
    InitializeGridFunctions();
    InterpolateDomainParameters();
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

void Domain_Parameters::InterpolateDomainParameters() {

    // interpolate domain parameter from distance function
    for (int vi = 0; vi < nV; vi++) {
        (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (Constants::zeta * Constants::dh))); // 
        // tanh(x) transitions -1 (electrolyte) to 1 (particle), adding 1.0 and multiplying 0.5 shifts range to 0 to 1
        // psi close to 1 is particle, psi close to 0 is electrolyte
        // large dsF = close to 1, so inside particle ; small or negative dsF = close to -1, so inside electrolyte

        (*pse)(vi) = 1.0 - (*psi)(vi);
        // pse close to 0 is particle, pse close to 1 is electrolyte

        (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (Constants::zeta * Constants::dh)), 2) - 1.0) / (2 * Constants::zeta * Constants::dh);
        // AvP is the rate of change of psi; used in reaction rates

        if ((*psi)(vi) < Constants::eps) { (*psi)(vi) = Constants::eps; }
        if ((*pse)(vi) < Constants::eps) { (*pse)(vi) = Constants::eps; }
    }

    AvB = std::make_unique<mfem::ParGridFunction>(*AvP);

    for (int vi = 0; vi < nV; vi++) {
        if ((*AvP)(vi) * Constants::dh < 1.0e-3) { (*AvP)(vi) = 0.0; }
        if ((*AvB)(vi) * Constants::dh < 1.0e-6) { (*AvB)(vi) = 0.0; }
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
    
    // Calculate totals for Psi and Pse fields
    CalculateTotalPhaseField(*psi, tPsi, gtPsi);
    CalculateTotalPhaseField(*pse, tPse, gtPse);

    // Compute the target current using Psi
    CalculateTargetCurrent(tPsi);
}

void Domain_Parameters::CalculateTargetCurrent(double total_psi) {
    
    // Compute target current based on total Psi, rho, Cr, and constants
    trgI = total_psi * Constants::rho * (0.9 - 0.3) / (3600.0 / Constants::Cr);

    // Perform global MPI reduction to get the total target current
    MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void Domain_Parameters::PrintInfo() {

    cout << "Total Psi: " << gtPsi << endl;
    cout << "Total Pse: " << gtPse << endl;
    cout << "Target Current: " << gTrgI << endl;

}
