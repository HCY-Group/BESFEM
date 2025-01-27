#ifndef MESH_BASE_HPP
#define MESH_BASE_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>

class MeshBase {
protected:
    std::unique_ptr<mfem::Mesh> globalMesh;             // Global serial mesh
    std::unique_ptr<mfem::ParMesh> parallelMesh;        // Parallel mesh
    std::shared_ptr<mfem::FiniteElementSpace> feSpace;  // Serial finite element space
    std::shared_ptr<mfem::FiniteElementSpace> globalfespace; // Parallel finite element space
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; // Parallel finite element space
    mfem::Vector elementVolumes;                       // Element volumes
    mfem::Array<int> boundaryMarkers;                  // Boundary markers

public:
    MeshBase();
    virtual ~MeshBase();

    // Mesh initialization
    void InitializeGlobalMesh(const char* meshFile);
    void InitializeParallelMesh(MPI_Comm comm);

    // Finite element space setup
    void SetupFiniteElementSpace(int order);

    // Parallel finite element space setup
    void SetupParFiniteElementSpace(int order);

    // // Element volume calculation
    // void CalculateElementVolumes();

    // // Boundary condition setup
    // virtual void SetupBoundaryConditions() = 0;

    // // Output functions
    // void OutputToParaview(const std::string &fileName, const std::string &varName, mfem::GridFunction *gf);

    // Print Mesh Information
    void PrintMeshInfo();

    // Getters for derived classes
    mfem::ParMesh *GetParallelMesh() const { return parallelMesh.get(); }
    mfem::ParFiniteElementSpace *GetParFiniteElementSpace() const { return parfespace.get(); }

    int nV; ///< Number of vertices
    int nE; ///< Number of elements
    int nC; ///< Number of corners per element

};

#endif // MESH_BASE_HPP
