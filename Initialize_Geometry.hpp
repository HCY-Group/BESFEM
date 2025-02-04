#ifndef INITIALIZE_GEOMETRY_HPP
#define INITIALIZE_GEOMETRY_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace std;


class Initialize_Geometry {
private:
    std::vector<std::vector<std::vector<int>>> tiffData;
    bool tiffDataLoaded = false;

protected:

    mfem::Vector elementVolumes;                       // Element volumes
    mfem::Array<int> boundaryMarkers;                  // Boundary markers
    std::vector<std::vector<std::vector<int>>> data;


public:
    Initialize_Geometry();
    virtual ~Initialize_Geometry();

    void InitializeMesh(const char* meshFile, MPI_Comm comm, int order);

    // Mesh initialization
    void InitializeGlobalMesh(const char* meshFile);
    void InitializeParallelMesh(MPI_Comm comm);
    std::vector<std::vector<std::vector<int>>>ReadTiffFile(const char* meshFile);
    std::unique_ptr<mfem::Mesh>CreateGlobalMeshFromTiffData(const std::vector<std::vector<std::vector<int>>>& tiffData);

    // Finite element space setup
    void SetupFiniteElementSpace(int order);

    // Parallel finite element space setup
    void SetupParFiniteElementSpace(int order);

    // Assign global values
    void AssignGlobalValues(const char* mesh_file);

    // Map global values to local
    void MapGlobalToLocal(const char* meshFile);

    // // Boundary condition setup
    // virtual void SetupBoundaryConditions() = 0;

    // // Output functions
    // void OutputToParaview(const std::string &fileName, const std::string &varName, mfem::GridFunction *gf);

    // Print Mesh Information
    void PrintMeshInfo();

    // Getters for derived classes
    mfem::ParMesh *GetParallelMesh() const { return parallelMesh.get(); }
    std::shared_ptr<mfem::ParFiniteElementSpace> GetParFiniteElementSpace() const {
        return parfespace;
    }
    int nV; ///< Number of vertices
    int nE; ///< Number of elements
    int nC; ///< Number of corners per element

    int gei;                // global element indices
    int ei;                 // local element indices

    std::unique_ptr<mfem::Mesh> globalMesh;             // Global serial mesh
    std::unique_ptr<mfem::ParMesh> parallelMesh;        // Parallel mesh
    std::shared_ptr<mfem::FiniteElementSpace> feSpace;  // Serial finite element space
    std::shared_ptr<mfem::FiniteElementSpace> globalfespace; // Parallel finite element space
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; // Parallel finite element space

    double Onm; ///< Number of grid function entries
    std::unique_ptr<mfem::GridFunction> gDsF; ///< Global distance function grid
    std::unique_ptr<mfem::ParGridFunction> dsF; ///< distance function grid
    std::unique_ptr<mfem::GridFunction> gVox; ///< Global vox function grid
    std::unique_ptr<mfem::ParGridFunction> Vox; ///< Vox function grid


};

#endif // INITIALIZE_GEOMETRY_HPP
