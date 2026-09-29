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
 * - Reading TIFF voxel data and constructing Cartesian meshes
 * - Creating global and parallel MFEM meshes
 * - Initializing serial and parallel FE spaces (H1 and DG)
 * - Mapping global fields to parallel fields
 * - Filtering solid, electrolyte, and particle-group voxel masks
 * - Refining phase interfaces and updating parallel spaces
 *
 * Boundary-condition setup and potential anchoring are handled by
 * BoundaryConditions in the current simulation driver.
 *
 * It serves as the lowest-level mesh/geometry infrastructure on which all
 * physics components (CnA, CnC, CnE, potentials, reactions) depend.
 */
class Initialize_Geometry {
private:
    const SimulationConfig& cfg; ///< Borrowed configuration; must outlive this geometry handler.

    /**
     * @brief Refine the half-cell interface region using a temporary solid mask.
     * Uses cfg.amr_levels and updates parallel spaces after refinement.
     * Returns immediately when no AMR levels are requested.
     */
    void HalfCellAMR();
    /**
     * @brief Refine full-cell interfaces using temporary anode and cathode masks.
     * Uses cfg.amr_levels and updates parallel spaces after refinement.
     * Returns immediately when no AMR levels are requested.
     */
    void FullCellAMR();

    /**
     * @brief Update parallel H1, scalar DG, vector DG, and any existing Vox field.
     * @pre Parallel finite-element spaces have been initialized.
     */
    void UpdateSpacesAfterAMR();

    /**
     * @brief Refresh local vertex, element, and element-corner counts.
     * Sets nC to zero when this rank has no elements.
     */
    void UpdateMeshData();

    /**
     * @brief Allocate and zero solid, electrolyte, and per-group half-cell masks.
     * @pre parfespace and particle_labels are prepared.
     */
    void AllocateHalfCellGeometryFields();
    /**
     * @brief Build half-cell masks on the current parallel mesh.
     * Builds collector-connected particle masks or copies the total solid mask
     * when cfg.combine_particle_groups is enabled.
     * @pre Half-cell geometry fields have been allocated.
     */
    void BuildHalfCellGeometryFields();

    /**
     * @brief Allocate and zero total and per-group full-cell phase masks.
     * @pre parfespace and both electrode label lists are prepared.
     */
    void AllocateFullCellGeometryFields();
    /**
     * @brief Filter full-cell anode, cathode, electrolyte, and particle masks.
     * Per-label calls do not request boundary-connectivity pruning.
     * @pre Full-cell geometry fields have been allocated.
     */
    void BuildFullCellGeometryFields();

    /**
     * @brief Reduce mesh statistics across ranks and print them on rank zero.
     * @param level Refinement level displayed with element counts and size bounds.
     * @note All MPI ranks must participate.
     */
    void PrintAMRMeshInfo(int level) const;

    /**
     * @brief Smooth a binary voxel mask and project it into the output space.
     *
     * Projects mask values as -1/+1 into DG, applies a PDE filter with weight
     * 3 * cfg.dh, rescales by (value + 1) / 2, and synchronizes true DOFs.
     * @param mask Flattened binary mask, indexed as x + nx * (y + ny * z).
     * @param nx Number of voxel columns.
     * @param ny Number of voxel rows.
     * @param nz Number of slices; one for a 2D image.
     * @param[out] filt_gf Filtered mask on the parallel finite-element space.
     * @pre Mesh and parallel spaces exist; mask contains nx * ny * nz entries.
     */
    void ApplyPDEFilterToMask(const std::vector<uint8_t>& mask, int nx, int ny, int nz, mfem::ParGridFunction& filt_gf);


protected:
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
     * @brief Initialize meshes, spaces, and filtered masks for a half cell.
     *
     * Reads TIFF data, creates serial and parallel spaces, maps voxel labels,
     * identifies particle groups, applies AMR, and builds final phase masks.
     * Uses cfg.half_electrode to select the electrode and writes geometry files.
     *
     * @param meshFile Path to a TIFF file with the .tif extension.
     * @param comm MPI communicator.
     * @param order Polynomial order for FE space.
     */
    void InitializeMesh(const char* meshFile, MPI_Comm comm, int order);


    /**
     * @brief Initialize meshes, spaces, and filtered masks for a full cell.
     *
     * Merges signed electrode TIFF stacks, creates meshes and spaces, maps
     * labels, applies AMR, and builds total and per-group phase masks.
     * Writes the merged geometry preview and parallel geometry fields.
     *
     * @param AnodeMeshFile TIFF stack with nonpositive anode labels.
     * @param CathodeMeshFile TIFF stack with nonnegative cathode labels.
     * @param comm MPI communicator.
     * @param order Polynomial order for FEspace.
     */
    void InitializeMesh(const char* AnodeMeshFile, const char* CathodeMeshFile, MPI_Comm comm, int order);

    /**
     * @brief Concatenate signed anode and cathode TIFF stacks along x.
     *
     * @param AnodeMeshFile TIFF stack with labels <= 0.
     * @param CathodeMeshFile TIFF stack with labels >= 0.
     * @return Merged voxel array indexed by [z][y][x], without added separator columns.
     * @throws std::runtime_error If either stack is empty, row/depth counts differ,
     *         or electrode labels have the wrong sign.
     * @note Writes a PGM preview of the first merged slice.
     */
    std::vector<std::vector<std::vector<int>>> MergeMeshes(const char *AnodeMeshFile, const char *CathodeMeshFile);

    /**
     * @brief Sample TIFF labels onto the serial voxel field gVox.
     * Uses cfg.coarsen_factor to select voxel samples in Cartesian vertex order.
     * @pre Nonempty tiffData and a compatible serial space are available.
     * @throws std::runtime_error If globalfespace is not initialized.
     */
    void AssignGlobalValues();
    /**
     * @brief Populate Vox from gVox using local-to-global element indices.
     * Also initializes mesh counts and vertex-index workspaces.
     * @pre Serial and parallel meshes, spaces, and gVox are initialized.
     * @throws std::runtime_error If a required mesh is missing.
     */
    void MapGlobalToLocal();
    /**
     * @brief Extract signed electrode labels from tiffData.
     * Negative labels belong to the anode and are ordered by absolute value;
     * positive cathode labels are ascending. Zero denotes electrolyte.
     * When combining groups, replaces nonempty lists with -1 and +1 respectively.
     */
    void FullCellParticleLabels();


    /**
     * @brief Initialize the global mesh from voxel data.
     * @param voxelData Nonempty voxel array indexed by [z][y][x].
     * @throws std::invalid_argument If voxel data is empty.
     * @throws std::runtime_error If the generated mesh has no elements.
     * @note Copies voxelData into tiffData and enables nonconforming refinement.
     */
    void InitializeGlobalMesh(const std::vector<std::vector<std::vector<int>>> &voxelData);

    /**
     * @brief Read a .tif file and construct a serial mesh supporting refinement.
     *
     * @param meshFile Path to a TIFF file with the .tif extension.
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
     * @return Cropped integer array indexed by [page][row][column].
     * @note Uses crop bounds and particle-color settings from cfg.
     */
    std::vector<std::vector<std::vector<int>>> ReadTiffFile(const char* meshFile);

    /**
     * @brief Construct global mesh from voxelized TIFF data.
     *
     * @param tiffData Nonempty rectangular voxel array indexed by [z][y][x].
     * @pre cfg.coarsen_factor is positive and produces at least one element
     *      in each active dimension.
     * @note One slice creates quadrilaterals; multiple slices create hexahedra.
     *       Element counts use cfg.coarsen_factor; physical extents use cfg.dh.
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
     * @brief Report when the parallel mesh has not been initialized.
     * @note Currently prints nothing when the parallel mesh exists.
     */
    void PrintMeshInfo();


    /**
     * @brief Save the first voxel slice as an 8-bit binary PGM preview.
     *
     * @param data Voxel array indexed by [z][y][x]; only data[0] is saved.
     * @param filename Complete output filename.
     * @note Maps the slice label range to 0-255; a constant slice becomes zero.
     *       Reports empty input or file-open failure and returns.
     */
    void SaveTiffDataToPGM(const std::vector<std::vector<std::vector<int>>> &data,
                       const std::string &filename);

    /**
     * @brief Build a connectivity-pruned phase mask and apply the PDE filter.
     *
     * Solid is positive in half cells, negative for full-cell anodes, and
     * positive for full-cell cathodes. Zero denotes electrolyte. Retains solid
     * connected to its collector, half-cell electrolyte connected to the
     * opposite boundary, or full-cell electrolyte touching both electrodes.
     * @note All MPI ranks participate in mask broadcast and filtering.
     * @param[out] filt_gf Filtered phase mask.
     * @param phase Geometry phase (SOLID or ELECTROLYTE).
     * @param cell_mode Cell mode (HALF or FULL).
     * @param electrode ANODE or CATHODE; selects the solid phase and collector side.
     */
    void ComputePDEFilter(
        mfem::ParGridFunction &filt_gf,
        sim::GeometryPhase phase,
        sim::CellMode cell_mode,
        sim::Electrode electrode);

    /**
     * @brief Compute a filtered mask for a particle label or combined group.
     *
     * Builds the electrode solid network, optionally retains its boundary-connected
     * component, then selects the requested label and filters the mask. When
     * cfg.combine_particle_groups is enabled, selects all solids of the electrode
     * instead of matching target_label. All MPI ranks participate.
     *
     * @param[out] filt_gf Filtered particle-group mask.
     * @param target_label Voxel label to isolate.
     * @param keep_boundary_connected Whether to keep only the boundary-connected region.
     * @param seed_side_or_face Seed boundary used when keep_boundary_connected is true.
     * @param cell_mode HALF or FULL cell mode.
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
     * @return Sorted unique nonzero labels; zero (electrolyte) is excluded.
     */
    std::vector<int> GetParticleLabelsFromTiff() const;

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * @brief Access distributed mesh.
     * @return Borrowed parallel mesh pointer, or nullptr before initialization.
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
    int nV = 0; ///< Number of vertices on this MPI rank.
    int nE = 0; ///< Number of elements on this MPI rank.
    int nC = 0; ///< Corners per element.

    int gei = 0; ///< Global element index.
    int ei  = 0; ///< Local element index.


    mfem::Array<int> gVTX; ///< Global vertex IDs of current element.
    mfem::Array<int> VTX;  ///< Local vertex IDs of current element.

    std::unique_ptr<mfem::Mesh> globalMesh; ///< Global serial mesh.
    mfem::Array<HYPRE_BigInt> E_L2G;         ///< Local-to-global element mapping.


    int myid = 0; ///< MPI rank.

    std::shared_ptr<mfem::ParMesh> parallelMesh; ///< Distributed parallel mesh.

    std::shared_ptr<mfem::FiniteElementSpace> globalfespace; ///< Global serial finite element space.

    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; ///< Parallel H1 finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace_dg; ///< Parallel DG finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> pardimfespace_dg; ///< Vector-valued parallel DG finite element space.


    std::unique_ptr<mfem::GridFunction> gVox; ///< Global voxel-label field.
    std::unique_ptr<mfem::ParGridFunction> Vox; ///< Parallel voxel-label field.

    std::vector<std::vector<std::vector<int>>> tiffData; ///< TIFF voxel labels indexed by [z][y][x].

    std::unique_ptr<mfem::H1_FECollection> gfec; ///< Serial H1 finite element collection.
    std::unique_ptr<mfem::H1_FECollection> pfec; ///< Parallel H1 finite element collection.
    std::unique_ptr<mfem::DG_FECollection> pfec_dg; ///< Parallel DG finite element collection.

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
