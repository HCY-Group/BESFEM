#include "MeshBase.hpp"
#include "Constants.hpp"
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

MeshBase::MeshBase() 
    // : mesh_file(Constants::mesh_file), dsF_file(Constants::dsF_file), order(Constants::order), dh(Constants::dh), 
    //   eps(Constants::eps), rho(Constants::rho), Cr(Constants::Cr), ze(Constants::ze)


{}

MeshBase::~MeshBase() {}

void MeshBase::InitializeGlobalMesh(const char* meshFile) {
    globalMesh = std::make_unique<mfem::Mesh>(meshFile);
    globalMesh->EnsureNCMesh(true);
}

void MeshBase::InitializeParallelMesh(MPI_Comm comm) {
    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before creating a parallel mesh.");
    }
    parallelMesh = std::make_unique<mfem::ParMesh>(comm, *globalMesh);
}

void MeshBase::SetupFiniteElementSpace(int order) {
    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before setting up FE space.");
    }
    auto fec = new mfem::H1_FECollection(order, globalMesh->Dimension());
    globalfespace = std::make_shared<mfem::FiniteElementSpace>(globalMesh.get(), fec);
}

void MeshBase::SetupParFiniteElementSpace(int order) {
    if (!parallelMesh) {
        throw std::runtime_error("Parallel mesh must be initialized before setting up FE space.");
    }
    auto fec = new mfem::H1_FECollection(order, parallelMesh->Dimension());
    parfespace = std::make_shared<mfem::ParFiniteElementSpace>(parallelMesh.get(), fec);
}

// void MeshBase::CalculateElementVolumes() {
//     if (!parallelMesh) {
//         throw std::runtime_error("Parallel mesh must be initialized before calculating element volumes.");
//     }
//     int numElements = parallelMesh->GetNE();
//     elementVolumes.SetSize(numElements);
//     for (int i = 0; i < numElements; ++i) {
//         elementVolumes[i] = parallelMesh->GetElementVolume(i);
//     }
// }

// void MeshBase::MapGlobalToLocal()

// void MeshBase::OutputToParaview(const std::string &fileName, const std::string &varName, mfem::GridFunction *gf) {
//     auto pd = new mfem::ParaViewDataCollection(fileName, gf->FESpace()->GetMesh());
//     pd->RegisterField(varName, gf);
//     pd->SetLevelsOfDetail(1);
//     pd->SetDataFormat(mfem::VTKFormat::BINARY);
//     pd->SetHighOrderOutput(true);
//     pd->Save();
//     delete pd;
// }

void MeshBase::PrintMeshInfo() {
    if (!parallelMesh) {
        std::cout << "Parallel mesh not initialized.\n";
        return;
    }

    nV = parallelMesh->GetNV();					// number of vertices
	nE = parallelMesh->GetNE();					// number of elements
	nC = pow(2, parallelMesh->Dimension());		// number of corner vertices

    std::cout << "Number of vertices: " << nV << "\n";
    std::cout << "Number of elements: " << nE << "\n";
}
