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
    : geometry(geo), nV(geo.nV), nE(geo.nE), nC(geo.nC), dsF(std::move(geo.dsF)), pmesh(std::move(geo.parallelMesh))
{}

// Destructor
Domain_Parameters::~Domain_Parameters() {}

void Domain_Parameters::SetupDomainParameters(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace){
    
    InitializeGridFunctions(fespace);
    InterpolateDomainParameters(fespace);
    CalculatePhasePotentialsAndTargetCurrent();

    PrintInfo();

    // pmesh->Save("pmesh");
    // dsF->Save("dsF_DP");
    // psi->Save("psi");
    // pse->Save("pse");
}

void Domain_Parameters::InitializeGridFunctions(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace) {

    if (!fespace) {
        throw std::runtime_error("Finite element space is not initialized.");
    }
    
    psi = make_unique<mfem::ParGridFunction>(fespace.get());
    pse = make_unique<mfem::ParGridFunction>(fespace.get());
    AvP = make_unique<mfem::ParGridFunction>(fespace.get());
    AvB = make_unique<mfem::ParGridFunction>(fespace.get());


}

void Domain_Parameters::InterpolateDomainParameters(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace) {

    // interpolate domain parameter from distance function
    for (int vi = 0; vi < nV; vi++) {
        // (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (zeta * dh))); // must use * to dereference to access
        (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (Constants::ze))); // must use * to dereference to access

        (*pse)(vi) = 1.0 - (*psi)(vi);
        // (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (zeta * dh)), 2) - 1.0) / (2 * zeta * dh);
        (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (Constants::ze)), 2) - 1.0) / (2 * Constants::ze);


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
