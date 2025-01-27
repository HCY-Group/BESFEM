// This Mesh_Test code tests handling the mesh as done in Anna's OOP, but now using Mesh Base class MeshBase

#include "Mesh_Test.hpp"
#include "Constants.hpp"
#include <iostream>
#include <cmath>
#include <mpi.h>

Mesh_Test::Mesh_Test() = default;

Mesh_Test::~Mesh_Test() = default;

void Mesh_Test::Initialize() {
    
    // Initialize the global mesh
    MeshBase::InitializeGlobalMesh(Constants::mesh_file);

    // Initialize the parallel mesh
    MeshBase::InitializeParallelMesh(MPI_COMM_WORLD);

    // Set up the finite element space
    MeshBase::SetupFiniteElementSpace(Constants::order);

    // Set up the parallel finite element space
    MeshBase::SetupParFiniteElementSpace(Constants::order);

    MeshBase::AssignGlobalValues(Constants::mesh_file, data);

    MeshBase::MapGlobalToLocal(Constants::mesh_file);

    // // Initialize grid functions
    // psi = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    // pse = std::make_unique<mfem::ParGridFunction>(parFeSpace.get());

    // // Compute element volumes
    // CalculateElementVolumes();

    PrintMeshInfo();
}

// void Mesh_Test::SetupBoundaryConditions() {
//     if (!parallelMesh) {
//         throw std::runtime_error("Parallel mesh must be initialized before setting up boundary conditions.");
//     }

//     // Example: Setup Neumann BCs on the west boundary
//     boundaryMarkers.SetSize(parallelMesh->bdr_attributes.Max());
//     boundaryMarkers = 0;
//     boundaryMarkers[3] = 1; // Apply Neumann BC to west boundary
// }

// void Mesh_Test::CalculatePhaseFields() {
//     if (!psi || !pse) {
//         throw std::runtime_error("Grid functions must be initialized before calculating phase fields.");
//     }

//     double zeta = 0.1; // Example parameter
//     double dh = 0.01;  // Example parameter

//     for (int i = 0; i < psi->Size(); ++i) {
//         double distance = (*psi)(i); // Use psi for intermediate calculations
//         (*psi)(i) = 0.5 * (1.0 + tanh(distance / (zeta * dh)));
//         (*pse)(i) = 1.0 - (*psi)(i);
//     }

//     // Sum up totals across elements
//     for (int i = 0; i < elementVolumes.Size(); ++i) {
//         totalPsi += (*psi)(i) * elementVolumes(i);
//         totalPse += (*pse)(i) * elementVolumes(i);
//     }

//     // Perform global MPI reduction
//     MPI_Allreduce(MPI_IN_PLACE, &totalPsi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, &totalPse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//     // Compute target current
//     double rho = 1.0;   // Example density parameter
//     double Cr = 1.0;    // Example reaction constant
//     targetCurrent = totalPsi * rho * (0.9 - 0.3) / (3600.0 / Cr);
//     MPI_Allreduce(MPI_IN_PLACE, &targetCurrent, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// }

// void Mesh_Test::PrintResults() const {
//     std::cout << "Total Psi: " << totalPsi << "\n";
//     std::cout << "Total Pse: " << totalPse << "\n";
//     std::cout << "Target Current: " << targetCurrent << "\n";
// }
