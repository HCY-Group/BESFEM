#ifndef INITIALIZE_GEOMETRY_HPP
#define INITIALIZE_GEOMETRY_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>
#include "SimTypes.hpp"

using namespace std;

/**
 * @class Initialize_Geometry
 * @brief Central geometry handler for BESFEM simulations.
 *
 * This class is responsible for:
 * - Loading meshes (serial or voxel-based)
 * - Creating global and parallel MFEM meshes
 * - Initializing serial and parallel FE spaces (H1 and DG)
 * - Mapping global fields to parallel fields
 * - Loading distance functions ψ (solid), ψ_e (electrolyte), and voxel data
 * - Setting up boundary-condition markers
 * - Selecting pinned DOFs for potential anchoring
 *
 * It serves as the lowest-level mesh/geometry infrastructure on which all
 * physics components (CnA, CnC, CnE, potentials, reactions) depend.
 */
class Initialize_Geometry {
private:
    bool tiffDataLoaded = false; ///< Whether TIFF voxel data has been loaded.

protected:
    mfem::Vector elementVolumes;     ///< Per-element volumes (global or parallel).
    mfem::Array<int> boundaryMarkers; ///< Boundary attribute markers.
    std::vector<std::vector<std::vector<int>>> data; ///< Raw voxel data container.

public:
    Initialize_Geometry();
    virtual ~Initialize_Geometry();

    // -------------------------------------------------------------------------
    // Distance function preprocessing
    // -------------------------------------------------------------------------

    /**
     * @brief Adjust distance function file to account for mesh spacing (dh).
     *
     * Applies smoothing, scaling, or clipping to ensure consistent SBM interface
     * thickness across meshes of different resolutions.
     *
     * @param distanceFile Path to the distance function file.
     * @param mesh_type Type of mesh ("ml" = MATLAB, "v" = voxel).
     */
    void AdjustDistanceFile(const char* distanceFile, const char* mesh_type);

    // -------------------------------------------------------------------------
    // Mesh initialization (serial + parallel)
    // -------------------------------------------------------------------------

    /**
     * @brief Initialize mesh and distance fields for a half-cell simulation.
     *
     * Creates the global mesh, assigns the distance field ψ or ψ_e from file,
     * builds FE spaces, and initializes parallel mesh/FESpaces.
     *
     * @param meshFile Path to mesh file.
     * @param distanceFile Distance function file.
     * @param mesh_type Type of mesh ("ml" = MATLAB, "v" = voxel).
     * @param comm MPI communicator.
     * @param order Polynomial order for FE space.
     */
    void InitializeMesh(const char* meshFile, const char* distanceFile, const char* mesh_type, MPI_Comm comm, int order);

    /**
     * @brief Initialize mesh and distance fields for a full-cell simulation.
     *
     * Loads separate anode and cathode distance functions.
     *
     * @param meshFile Path to mesh file.
     * @param distanceFileA Anode distance function file (ψ_A).
     * @param distanceFileC Cathode distance function file (ψ_C).
     * @param mesh_type Type of mesh ("ml" = MATLAB, "v" = voxel).
     * @param comm MPI communicator.
     * @param order Polynomial order for FE space.
     */
    void InitializeMesh(const char* meshFile, const char* distanceFileA, const char* distanceFileC, const char* mesh_type, MPI_Comm comm, int order);

    /**
     * @brief Load and construct global serial MFEM mesh.
     *
     * @param meshFile Path to mesh file.
     */
    void InitializeGlobalMesh(const char* meshFile);

    /**
     * @brief Build a distributed (parallel) mesh from the global mesh.
     *
     * @param comm MPI communicator.
     */
    void InitializeParallelMesh(MPI_Comm comm);

    // -------------------------------------------------------------------------
    // Voxel / TIFF mesh support
    // -------------------------------------------------------------------------

    /**
     * @brief Read TIFF voxel data for voxel-mesh construction.
     *
     * @param meshFile TIFF volume file.
     * @return 3D voxel array.
     */
    std::vector<std::vector<std::vector<int>>> ReadTiffFile(const char* meshFile);

    /**
     * @brief Construct global mesh from voxelized TIFF data.
     *
     * @param tiffData 3D voxel grid.
     * @return Newly constructed MFEM mesh.
     */
    std::unique_ptr<mfem::Mesh>
    CreateGlobalMeshFromTiffData(const std::vector<std::vector<std::vector<int>>>& tiffData);

    // -------------------------------------------------------------------------
    // Finite element spaces
    // -------------------------------------------------------------------------

    /**
     * @brief Setup global serial FE space.
     *
     * @param order FE polynomial order.
     */
    void SetupFiniteElementSpace(int order);

    /**
     * @brief Setup parallel FE space.
     *
     * @param order FE polynomial order.
     */
    void SetupParFiniteElementSpace(int order);

    // -------------------------------------------------------------------------
    // Global → parallel mapping and distance functions
    // -------------------------------------------------------------------------

    /**
     * @brief Assign global fields: distance functions, voxel mask, etc.
     *
     * Loads, projects, and assigns global fields onto the global mesh.
     *
     * @param mesh_file Mesh path.
     * @param distanceFile Distance field file.
     * @param gDsF_out Output global distance function ψ.
     */
    void AssignGlobalValues(const char* mesh_file, const char* distanceFile, std::unique_ptr<mfem::GridFunction>& gDsF_out);

    /**
     * @brief Transfer global mesh quantities to parallel fields.
     *
     * @param meshFile Mesh file for consistency checking.
     */
    void MapGlobalToLocal(const char* meshFile);

    // -------------------------------------------------------------------------
    // Boundary conditions and pinning
    // -------------------------------------------------------------------------

    /**
     * @brief Set Neumann/Dirichlet boundary condition markers.
     *
     * Uses boundary attributes + cell mode to mark anode, cathode, electrolyte regions.
     *
     * @param mode HALF or FULL cell.
     * @param electrode ANODE, CATHODE, BOTH.
     */
    void SetupBoundaryConditions(sim::CellMode mode, sim::Electrode electrode);

    /**
     * @brief Print mesh information (elements, vertices, dimensions).
     */
    void PrintMeshInfo();

    /**
     * @brief Select and assign a pinned DOF for potential anchoring.
     *
     * @param fespace FE space in which to anchor a single true DOF.
     */
    void SetupPinnedDOF(mfem::ParFiniteElementSpace& fespace);

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * @brief Access distributed mesh.
     * @return Pointer to parallel mesh.
     */
    mfem::ParMesh *GetParallelMesh() const { return parallelMesh.get(); }

    /**
     * @brief Access parallel FE space.
     * @return Parallel FE space (shared_ptr).
     */
    std::shared_ptr<mfem::ParFiniteElementSpace>
    GetParFiniteElementSpace() const { return parfespace; }

    // -------------------------------------------------------------------------
    // Public geometry/mesh fields
    // -------------------------------------------------------------------------
    int nV = 0; ///< Number of vertices.
    int nE = 0; ///< Number of elements.
    int nC = 0; ///< Corners per element.

    int gei = 0; ///< Global element index.
    int ei  = 0; ///< Local element index.

    // Boundary condition marker arrays
    mfem::Array<int> nbc_w_bdr, nbc_s_bdr, nbc_e_bdr, nbc_n_bdr;
    mfem::Array<int> nbc_bdr, dbc_bdr;
    mfem::Array<int> dbc_w_bdr, dbc_e_bdr;

    mfem::Array<int> ess_tdof_list_w, ess_tdof_list_e;

    mfem::Array<int> gVTX; ///< Global vertex IDs of current element.
    mfem::Array<int> VTX;  ///< Local vertex IDs of current element.

    // Mesh + FE space pointers
    std::unique_ptr<mfem::Mesh> globalMesh; ///< Serial/global mesh.
    std::shared_ptr<mfem::ParMesh> parallelMesh;

    std::shared_ptr<mfem::FiniteElementSpace> feSpace;
    std::shared_ptr<mfem::FiniteElementSpace> globalfespace;

    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace;
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace_dg;
    std::shared_ptr<mfem::ParFiniteElementSpace> pardimfespace_dg;

    mfem::Array<HYPRE_BigInt> E_L2G; ///< Local-to-global element mapping.

    double Onm = 0.0; ///< Number of grid function entries.

    // Distance fields
    std::unique_ptr<mfem::GridFunction> gDsF;
    std::unique_ptr<mfem::ParGridFunction> dsF;

    std::unique_ptr<mfem::GridFunction> gDsF_A, gDsF_C;
    std::unique_ptr<mfem::ParGridFunction> dsF_A, dsF_C;

    std::unique_ptr<mfem::GridFunction> gVox;
    std::unique_ptr<mfem::ParGridFunction> Vox;

    // TIFF voxel storage
    std::vector<std::vector<std::vector<int>>> tiffData;

    // FE collections
    std::unique_ptr<mfem::H1_FECollection> gfec, pfec;
    std::unique_ptr<mfem::DG_FECollection> pfec_dg;

    // Pinned DOF information
    mfem::Array<int> ess_tdof_potE;
    bool anchor_set = false;
    HYPRE_BigInt global_anchor_potE = -1;
    int anchor_owner_potE = -1;
    bool pin = false;

    mfem::Array<int> ess_tdof_listPinned;
    mfem::Array<int> boundary_dofs;
    mfem::Array<int> ess_tdof_marker;

    int myid = 0; ///< MPI rank.
    int rkpp = -1; ///< Rank that owns the pinned DOF.
};

#endif // INITIALIZE_GEOMETRY_HPP
