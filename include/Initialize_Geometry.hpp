#ifndef INITIALIZE_GEOMETRY_HPP
#define INITIALIZE_GEOMETRY_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>
#include "SimTypes.hpp"
#include "SimulationConfig.hpp"
#include <set>

/**
 * @file Initialize_Geometry.hpp
 * @brief Defines mesh, geometry, voxel, distance-field, and FE-space setup for BESFEM.
 */

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
    const SimulationConfig& cfg;

    void HalfCellAMR();
    void FullCellAMR();

    void UpdateSpacesAfterAMR();

    void UpdateMeshData();

    void AllocateHalfCellGeometryFields();
    void BuildHalfCellGeometryFields();

    void AllocateFullCellGeometryFields();
    void BuildFullCellGeometryFields();

    void PrintAMRMeshInfo(int level) const;

    void ApplyPDEFilterToMask(const std::vector<uint8_t>& mask, int nx, int ny, int nz, mfem::ParGridFunction& filt_gf);


protected:
    mfem::Vector elementVolumes;     ///< Per-element volumes (global or parallel).
    mfem::Array<int> boundaryMarkers; ///< Boundary attribute markers.
    std::vector<std::vector<std::vector<int>>> data; ///< Raw voxel data container.

public:

    /**
     * @brief Construct the geometry handler.
     *
     * Stores the simulation configuration and prepares geometry-related storage.
     *
     * @param cfg Reference to the simulation configuration.
     */
    Initialize_Geometry(const SimulationConfig& cfg);

    /// Destructor.
    virtual ~Initialize_Geometry();

    bool combine_particle_groups = false; ///< Whether to combine particle groups for performance.

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
     * @param comm MPI communicator.
     * @param order Polynomial order for FE space.
     */
    void InitializeMesh(const char* meshFile, MPI_Comm comm, int order);


    /**
     * @brief Initialize mesh and distance fields for a full-cell simulation.
     *
     * Loads separate anode and cathode distance functions.
     *
     * @param AnodeMeshFile Path to anode mesh file.
     * @param CathodeMeshFile Path to cathode mesh file.
     * @param comm MPI communicator.
     * @param order Polynomial order for FEspace.
     */
    void InitializeMesh(const char* AnodeMeshFile, const char* CathodeMeshFile, MPI_Comm comm, int order);

    /*
     * @brief Merge anode and cathode meshes.
     *
     * @param anode_mesh_file Path to anode mesh file.
     * @param cathode_mesh_file Path to cathode mesh file.
     * @return Merged mesh.
     */
    std::vector<std::vector<std::vector<int>>> MergeMeshes(const char *AnodeMeshFile, const char *CathodeMeshFile);

    void AssignGlobalValues();
    void MapGlobalToLocal();
    void FullCellParticleLabels();


    /**
     * @brief Initialize the global mesh from voxel data.
     * @param voxelData vector of voxel data.
     */
    void InitializeGlobalMesh(const std::vector<std::vector<std::vector<int>>> &voxelData);

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
    
    /**
     * @brief Save TIFF voxel data to a series of PGM files.
     * 
     * @param data voxel data.
     * @param filename Base filename for output PGM files.
     */
    void SaveTiffDataToPGM(const std::vector<std::vector<std::vector<int>>> &data,
                       const std::string &filename);

    /**
     * @brief Use MFEM-based solvers to compute distance function from voxel mask.
     * @param filt_gf Output filtered level set function.
     * @param phase Geometry phase (SOLID or ELECTROLYTE).
     * @param cell_mode Cell mode (HALF or FULL).
     * @param electrode Electrode type (ANODE, CATHODE, BOTH).
     */
    void ComputePDEFilter(
        mfem::ParGridFunction &filt_gf,
        sim::GeometryPhase phase,
        sim::CellMode cell_mode,
        sim::Electrode electrode);

    /**
     * @brief Compute a filtered distance field for a specific voxel label.
     *
     * Builds a label-specific mask from the TIFF data, applies the PDE filter,
     * and computes a corresponding distance field. This is used for individual
     * particle or material-region masks in segmented microstructures.
     *
     * @param filt_gf Output filtered level-set field.
     * @param target_label Voxel label to isolate.
     * @param keep_boundary_connected Whether to keep only the boundary-connected region.
     * @param seed_side_or_face Optional boundary side/face used as the seed region.
     * @param cell_mode The mode of the cell (e.g., 2D or 3D).
     * @param electrode The electrode configuration.
     */
    void ComputePDEFilterLabel(mfem::ParGridFunction &filt_gf,
            int target_label,
            bool keep_boundary_connected,
            sim::BoundarySide seed_side_or_face,
            sim::CellMode cell_mode,
            sim::Electrode electrode);

    /**
     * @brief Return the unique particle/material labels found in the TIFF data.
     *
     * @return Vector of integer labels present in the voxelized geometry.
     */
    std::vector<int> GetParticleLabelsFromTiff() const;

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
    mfem::Array<int> nbc_w_bdr, nbc_s_bdr, nbc_e_bdr, nbc_n_bdr; ///< Neumann boundary markers for west, south, east, and north boundaries.
    mfem::Array<int> nbc_bdr, dbc_bdr;                           ///< Global Neumann and Dirichlet boundary markers.
    mfem::Array<int> dbc_w_bdr, dbc_e_bdr;                      ///< Dirichlet markers for west and east boundaries.

    mfem::Array<int> ess_tdof_list_w, ess_tdof_list_e;          ///< Essential true DOFs on west and east boundaries.

    mfem::Array<int> gVTX; ///< Global vertex IDs of current element.
    mfem::Array<int> VTX;  ///< Local vertex IDs of current element.

    std::unique_ptr<mfem::Mesh> globalMesh; ///< Global serial mesh.
    mfem::Array<HYPRE_BigInt> E_L2G;         ///< Local-to-global element mapping.

    double Onm = 0.0; ///< Number of grid function entries.

    // Pinned DOF information
    mfem::Array<int> ess_tdof_potE;          ///< Essential true DOFs used to anchor the electrolyte potential.
    bool anchor_set = false;                 ///< Whether the global anchor has been selected.
    HYPRE_BigInt global_anchor_potE = -1;    ///< Global ID of the pinned electrolyte potential DOF.
    int anchor_owner_potE = -1;              ///< MPI rank that owns the pinned DOF.
    bool pin = false;                        ///< Whether this MPI rank owns the pinned DOF.

    mfem::Array<int> ess_tdof_listPinned;    ///< Essential true DOF list for the pinned node.
    mfem::Array<int> boundary_dofs;          ///< Boundary true DOFs.
    mfem::Array<int> ess_tdof_marker;        ///< Marker array identifying essential DOFs.

    int myid = 0; ///< MPI rank.
    int rkpp = -1; ///< Rank that owns the pinned DOF.

    std::shared_ptr<mfem::ParMesh> parallelMesh; ///< Distributed parallel mesh.

    std::shared_ptr<mfem::FiniteElementSpace> feSpace; ///< Serial H1 finite element space.
    std::shared_ptr<mfem::FiniteElementSpace> globalfespace; ///< Global serial finite element space.

    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; ///< Parallel H1 finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace_dg; ///< Parallel DG finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> pardimfespace_dg; ///< Vector-valued parallel DG finite element space.

    std::unique_ptr<mfem::GridFunction> gDsF; ///< Global serial distance field.
    std::unique_ptr<mfem::ParGridFunction> dsF; ///< Parallel distance field.

    std::unique_ptr<mfem::GridFunction> gDsF_A, gDsF_C; ///< Global anode/cathode distance fields.
    std::unique_ptr<mfem::ParGridFunction> dsF_A, dsF_C; ///< Parallel anode/cathode distance fields.

    std::unique_ptr<mfem::GridFunction> gVox; ///< Global voxel-label field.
    std::unique_ptr<mfem::ParGridFunction> Vox; ///< Parallel voxel-label field.

    std::vector<std::vector<std::vector<int>>> tiffData; ///< Raw TIFF voxel labels.

    std::unique_ptr<mfem::H1_FECollection> gfec, pfec; ///< Serial/parallel H1 finite element collections.
    std::unique_ptr<mfem::DG_FECollection> pfec_dg; ///< Parallel DG finite element collection.

    std::unique_ptr<mfem::ParGridFunction> distMask; ///< Unsigned distance-to-mask field.
    std::unique_ptr<mfem::ParGridFunction> distMaskSigned; ///< Optional signed distance-to-mask field.
    std::unique_ptr<mfem::ParGridFunction> MaskFilter; ///< Filtered solid-mask level-set field.
    std::unique_ptr<mfem::ParGridFunction> MaskFilterPse; ///< Filtered electrolyte-mask level-set field.

    std::unique_ptr<mfem::ParGridFunction> MaskFilterAnode; ///< Filtered anode mask field.
    std::unique_ptr<mfem::ParGridFunction> MaskFilterCathode; ///< Filtered cathode mask field.

    std::vector<int> anode_particle_labels; ///< Unique anode particle labels.
    std::vector<int> cathode_particle_labels; ///< Unique cathode particle labels.

    std::vector<std::unique_ptr<mfem::ParGridFunction>>MaskFiltersAnode; ///< Per-label filtered anode mask fields.
    std::vector<std::unique_ptr<mfem::ParGridFunction>>MaskFiltersCathode; ///< Per-label filtered cathode mask fields.

    std::vector<int> particle_labels; ///< Unique particle/material labels from the TIFF geometry.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> MaskFilters; ///< Per-label filtered mask fields.

};

#endif // INITIALIZE_GEOMETRY_HPP
