#ifndef INITIALIZE_GEOMETRY_HPP
#define INITIALIZE_GEOMETRY_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace std;

/**
 * @class Initialize_Geometry
 * @brief Handles mesh, finite element space, and boundary condition initialization 
 *        for parallel finite element battery simulations.
 *
 * This class provides methods to load mesh data, set up serial and parallel
 * finite element spaces, assign global values, and configure boundary conditions.
 */
class Initialize_Geometry {
private:
    bool tiffDataLoaded = false; ///< Flag indicating whether TIFF mesh data has been loaded.

protected:

    mfem::Vector elementVolumes;                       ///< Volumes of mesh elements
    mfem::Array<int> boundaryMarkers;                  ///< Boundary markers for applying conditions
    std::vector<std::vector<std::vector<int>>> data;  ///< 3D voxel data container


public:
    Initialize_Geometry();
    virtual ~Initialize_Geometry();

    /**
     * @brief Adjusts the distance file to ensure scaling for dh.
     * @param distanceFile Path to the distance function file.
     */
    void AdjustDistanceFile(const char* distanceFile);

    /**
     * @brief Initializes mesh and associated distance functions for half cell.
     * @param meshFile Path to mesh file.
     * @param distanceFile Path to distance function file.
     * @param comm MPI communicator.
     * @param order Finite element polynomial order.
     */
    void InitializeMesh(const char* meshFile, const char* distanceFile, MPI_Comm comm, int order);

    /**
     * @brief Initializes mesh and associated distance functions for full cell.
     * @param meshFile Path to mesh file.
     * @param distanceFileA Path to distance function file for the anode.
     * @param distanceFileC Path to distance function file for the cathode.
     * @param comm MPI communicator.
     * @param order Finite element polynomial order.
     */
    void InitializeMesh(const char* meshFile, const char* distanceFileA, const char* distanceFileC, MPI_Comm comm, int order);

    /**
     * @brief Initializes a global serial mesh from a file.
     * @param meshFile Path to mesh file.
     */
    void InitializeGlobalMesh(const char* meshFile);

    /**
     * @brief Initializes a parallel mesh from the global mesh.
     * @param comm MPI communicator.
     */
    void InitializeParallelMesh(MPI_Comm comm);

    /**
     * @brief Reads TIFF image data for voxelized mesh creation.
     * @param meshFile Path to TIFF file.
     * @return 3D vector of voxel data.
     */
    std::vector<std::vector<std::vector<int>>>ReadTiffFile(const char* meshFile);

    /**
     * @brief Creates a global mesh from voxelized TIFF data.
     * @param tiffData 3D voxel data.
     * @return Unique pointer to a new mfem::Mesh.
     */
    std::unique_ptr<mfem::Mesh>CreateGlobalMeshFromTiffData(const std::vector<std::vector<std::vector<int>>>& tiffData);

    /**
     * @brief Sets up the serial finite element space.
     * @param order Polynomial order.
     */
    void SetupFiniteElementSpace(int order);

    /**
     * @brief Sets up the parallel finite element space.
     * @param order Polynomial order.
     */
    void SetupParFiniteElementSpace(int order);

    /**
     * @brief Assigns global values such as distance functions and voxel data.
     * @param mesh_file Path to mesh file.
     * @param distanceFile Path to distance function file.
     * @param gDsF_out Output unique pointer to the assigned global distance function.
     */
    void AssignGlobalValues(const char* mesh_file, const char* distanceFile, std::unique_ptr<mfem::GridFunction>& gDsF_out);

    /**
     * @brief Maps global values to local parallel data structures.
     * @param meshFile Path to mesh file.
     */
    void MapGlobalToLocal(const char* meshFile);

    /**
     * @brief Sets up boundary condition markers for Dirichlet and Neumann boundaries.
     */
    void SetupBoundaryConditions();

    /**
     * @brief Prints information about the mesh (elements, vertices, etc.).
     */
    void PrintMeshInfo();

    /**
     * @brief Returns the parallel mesh pointer.
     * @return Pointer to parallel mesh.
     */
    mfem::ParMesh *GetParallelMesh() const { return parallelMesh.get(); }

    /**
     * @brief Returns the parallel finite element space.
     * @return Shared pointer to mfem::ParFiniteElementSpace.
     */
    std::shared_ptr<mfem::ParFiniteElementSpace> GetParFiniteElementSpace() const {
        return parfespace;
    }

    int nV; ///< Number of vertices
    int nE; ///< Number of elements
    int nC; ///< Number of corners per element

    int gei; ///< Global element index.
    int ei;  ///< Local element index.

    mfem::Array<int> nbc_w_bdr; ///< West Neumann Boundary Conditions
    mfem::Array<int> nbc_s_bdr; ///< West Neumann Boundary Conditions
    mfem::Array<int> nbc_e_bdr; ///< West Neumann Boundary Conditions
    mfem::Array<int> nbc_n_bdr; ///< West Neumann Boundary Conditions

    
    mfem::Array<int> ess_tdof_list_w; ///< Total DOF West
    mfem::Array<int> ess_tdof_list_e; ///< Total DOF East

    mfem::Array<int> dbc_w_bdr; ///< West Dirichlet Boundary Conditions
    mfem::Array<int> dbc_e_bdr; ///< East Dirichlet Boundary Conditions
    
    std::unique_ptr<mfem::Mesh> globalMesh;              ///< Global serial mesh.
    std::unique_ptr<mfem::ParMesh> parallelMesh;         ///< Parallel mesh.
    std::shared_ptr<mfem::FiniteElementSpace> feSpace;   ///< Serial finite element space.
    std::shared_ptr<mfem::FiniteElementSpace> globalfespace; ///< Global finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; ///< Parallel finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace_dg; ///< Parallel DG finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> pardimfespace_dg; ///< Parallel DG FE space (dim).


    double Onm; ///< Number of grid function entries.
    std::unique_ptr<mfem::GridFunction> gDsF; ///< Global distance function.
    std::unique_ptr<mfem::ParGridFunction> dsF; ///< Parallel distance function.

    std::unique_ptr<mfem::GridFunction>       gDsF_A;  ///< Global distance function for anode.
    std::unique_ptr<mfem::GridFunction>       gDsF_C;  ///< Global distance function for cathode.
    std::unique_ptr<mfem::ParGridFunction>    dsF_A;    /// Parallel distance function for anode.
    std::unique_ptr<mfem::ParGridFunction>    dsF_C;    ///< Parallel distance function for cathode.
    
    std::unique_ptr<mfem::GridFunction> gVox; ///< Global voxel function.
    std::unique_ptr<mfem::ParGridFunction> Vox; ///< Parallel voxel function.

    std::vector<std::vector<std::vector<int>>> tiffData; ///< TIFF voxel data.
    std::unique_ptr<mfem::H1_FECollection> gfec; ///< Global H1 finite element collection.
    std::unique_ptr<mfem::H1_FECollection> pfec; ///< Parallel H1 finite element collection.
    std::unique_ptr<mfem::DG_FECollection> pfec_dg; ///< Parallel DG finite element collection.



};

#endif // INITIALIZE_GEOMETRY_HPP
