#include "../include/Initialize_Geometry.hpp"
#include "../include/Constants.hpp"
#include "../include/readtiff.h"
#include "../include/SimTypes.hpp"
#include "../include/dist_solver.hpp"

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace std;
using sim::CellMode;
using sim::Electrode;
using sim::GeometryPhase;
using sim::BoundarySide;


// Constructor
Initialize_Geometry::Initialize_Geometry(const SimulationConfig& cfg)
    : cfg(cfg)
{}

// Destructor
Initialize_Geometry::~Initialize_Geometry() {}

static void KeepOnlyConnectedToBoundary_2D(std::vector<uint8_t> &solid, int nx, int ny, bool eight_conn, bool seed_all_boundaries, BoundarySide seed_side)
{
    // seed_side: -1 = use all boundaries; 0=left, 1=right, 2=bottom, 3=top
    auto id = [nx](int i, int j){ return i + nx*j; };

    std::vector<uint8_t> keep(nx*ny, 0);
    std::queue<std::pair<int,int>> q;

    auto push = [&](int i, int j){
        if (i < 0 || i >= nx || j < 0 || j >= ny) return;
        int k = id(i,j);
        if (!solid[k] || keep[k]) return;
        keep[k] = 1;
        q.push({i,j});
    };

    // seeds
    if (seed_all_boundaries)
    {
        for (int i=0;i<nx;i++){ push(i,0); push(i,ny-1); }
        for (int j=0;j<ny;j++){ push(0,j); push(nx-1,j); }
    }
    else
    {
        if (seed_side == BoundarySide::WEST) for (int j=0;j<ny;j++) push(0,j);         // left
        if (seed_side == BoundarySide::EAST) for (int j=0;j<ny;j++) push(nx-1,j);      // right
        if (seed_side == BoundarySide::SOUTH) for (int i=0;i<nx;i++) push(i,0);         // bottom
        if (seed_side == BoundarySide::NORTH) for (int i=0;i<nx;i++) push(i,ny-1);      // top
    }

    const int di4[4] = { 1,-1, 0, 0};
    const int dj4[4] = { 0, 0, 1,-1};
    const int di8[8] = { 1,-1, 0, 0, 1, 1,-1,-1};
    const int dj8[8] = { 0, 0, 1,-1, 1,-1, 1,-1};

    while (!q.empty())
    {
        auto [i,j] = q.front(); q.pop();
        if (!eight_conn)
            for (int t=0;t<4;t++) push(i+di4[t], j+dj4[t]);
        else
            for (int t=0;t<8;t++) push(i+di8[t], j+dj8[t]);
    }

    // remove islands
    for (int k=0;k<nx*ny;k++) if (solid[k] && !keep[k]) solid[k] = 0;
}

static void KeepOnlyConnectedToBoundary_3D(std::vector<uint8_t> &solid, int nx, int ny, int nz, bool twenty_six_conn, bool seed_all_boundaries, BoundarySide seed_face)
{
    // seed_face: -1 = all faces
    // 0=xmin, 1=xmax, 2=ymin, 3=ymax, 4=zmin, 5=zmax
    auto id = [=](int i,int j,int k){ return i + nx*j + nx*ny*k; };

    std::vector<uint8_t> keep(nx*ny*nz, 0);
    std::queue<std::tuple<int,int,int>> q;

    auto push = [&](int i,int j,int k){
        if (i<0||i>=nx||j<0||j>=ny||k<0||k>=nz) return;
        int idx = id(i,j,k);
        if (!solid[idx] || keep[idx]) return;
        keep[idx] = 1;
        q.push({i,j,k});
    };

    auto seed_face_fn = [&](BoundarySide face){
        if (face==BoundarySide::WEST) for (int k=0;k<nz;k++) for (int j=0;j<ny;j++) push(0,j,k);         // xmin
        if (face==BoundarySide::EAST) for (int k=0;k<nz;k++) for (int j=0;j<ny;j++) push(nx-1,j,k);      // xmax
        if (face==BoundarySide::SOUTH) for (int k=0;k<nz;k++) for (int i=0;i<nx;i++) push(i,0,k);         // ymin
        if (face==BoundarySide::NORTH) for (int k=0;k<nz;k++) for (int i=0;i<nx;i++) push(i,ny-1,k);      // ymax
        if (face==BoundarySide::BOTTOM) for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) push(i,j,0);         // zmin
        if (face==BoundarySide::TOP) for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) push(i,j,nz-1);      // zmax
    };

    if (seed_all_boundaries)
    {
        seed_face_fn(BoundarySide::WEST);
        seed_face_fn(BoundarySide::EAST);
        seed_face_fn(BoundarySide::SOUTH);
        seed_face_fn(BoundarySide::NORTH);
        seed_face_fn(BoundarySide::BOTTOM);
        seed_face_fn(BoundarySide::TOP);
    }
    else
    {
        seed_face_fn(seed_face);
    }

    // BFS neighbors
    if (!twenty_six_conn)
    {
        const int di[6]={ 1,-1, 0, 0, 0, 0};
        const int dj[6]={ 0, 0, 1,-1, 0, 0};
        const int dk[6]={ 0, 0, 0, 0, 1,-1};
        while(!q.empty())
        {
            auto [i,j,k]=q.front(); q.pop();
            for(int t=0;t<6;t++) push(i+di[t], j+dj[t], k+dk[t]);
        }
    }
    else
    {
        while(!q.empty())
        {
            auto [i,j,k]=q.front(); q.pop();
            for(int dk=-1; dk<=1; dk++)
            for(int dj=-1; dj<=1; dj++)
            for(int di=-1; di<=1; di++)
            {
                if (di==0 && dj==0 && dk==0) continue;
                push(i+di, j+dj, k+dk);
            }
        }
    }

    // remove islands
    for (int idx=0; idx<nx*ny*nz; idx++)
        if (solid[idx] && !keep[idx]) solid[idx] = 0;
}

static void KeepOnlyElectrolyteTouchingBothElectrodes_2D(std::vector<uint8_t> &electrolyte, const std::vector<std::vector<std::vector<int>>> &labels,
    int nx, int ny, bool eight_conn)
{
    MFEM_VERIFY(
        labels.size() == 1,
        "KeepOnlyElectrolyteTouchingBothElectrodes_2D requires 2D TIFF data.");

    auto id = [nx](int i, int j)
    {
        return i + nx * j;
    };

    std::vector<uint8_t> visited(nx * ny, 0);
    std::vector<uint8_t> keep(nx * ny, 0);

    const int di4[4] = {1, -1, 0, 0};
    const int dj4[4] = {0, 0, 1, -1};

    const int di8[8] = {1, -1, 0, 0, 1, 1, -1, -1};
    const int dj8[8] = {0, 0, 1, -1, 1, -1, 1, -1};

    const int neighbor_count = eight_conn ? 8 : 4;
    const int *di = eight_conn ? di8 : di4;
    const int *dj = eight_conn ? dj8 : dj4;

    int kept_components = 0;
    int removed_components = 0;

    for (int j0 = 0; j0 < ny; ++j0)
    {
        for (int i0 = 0; i0 < nx; ++i0)
        {
            const int start = id(i0, j0);

            if (!electrolyte[start] || visited[start])
            {
                continue;
            }

            std::queue<std::pair<int, int>> q;
            std::vector<int> component;

            bool touches_anode = false;
            bool touches_cathode = false;

            visited[start] = 1;
            q.push({i0, j0});

            while (!q.empty())
            {
                const auto [i, j] = q.front();
                q.pop();

                component.push_back(id(i, j));

                for (int n = 0; n < neighbor_count; ++n)
                {
                    const int ni = i + di[n];
                    const int nj = j + dj[n];

                    if (ni < 0 || ni >= nx ||
                        nj < 0 || nj >= ny)
                    {
                        continue;
                    }

                    const int neighbor_label = labels[0][nj][ni];
                    const int neighbor_id = id(ni, nj);

                    if (neighbor_label < 0)
                    {
                        touches_anode = true;
                    }
                    else if (neighbor_label > 0)
                    {
                        touches_cathode = true;
                    }
                    else if (electrolyte[neighbor_id] &&
                             !visited[neighbor_id])
                    {
                        visited[neighbor_id] = 1;
                        q.push({ni, nj});
                    }
                }
            }

            if (touches_anode && touches_cathode)
            {
                for (const int idx : component)
                {
                    keep[idx] = 1;
                }

                ++kept_components;
            }
            else
            {
                ++removed_components;
            }
        }
    }

    electrolyte.swap(keep);
    std::cout << "[Full Cell Connectivity] Electrolyte components kept: " << kept_components << ", removed: " << removed_components << "\n";
}

static void KeepOnlyElectrolyteTouchingBothElectrodes_3D(std::vector<uint8_t> &electrolyte, const std::vector<std::vector<std::vector<int>>> &labels, int nx, int ny, int nz, bool twenty_six_conn)
{
    auto id = [=](int i, int j, int k)
    {
        return i + nx * j + nx * ny * k;
    };

    std::vector<uint8_t> visited(nx * ny * nz, 0);
    std::vector<uint8_t> keep(nx * ny * nz, 0);

    int kept_components = 0;
    int removed_components = 0;

    for (int k0 = 0; k0 < nz; ++k0)
    {
        for (int j0 = 0; j0 < ny; ++j0)
        {
            for (int i0 = 0; i0 < nx; ++i0)
            {
                const int start = id(i0, j0, k0);

                if (!electrolyte[start] || visited[start])
                {
                    continue;
                }

                std::queue<std::tuple<int, int, int>> q;
                std::vector<int> component;

                bool touches_anode = false;
                bool touches_cathode = false;

                visited[start] = 1;
                q.push({i0, j0, k0});

                while (!q.empty())
                {
                    const auto [i, j, k] = q.front();
                    q.pop();

                    component.push_back(id(i, j, k));

                    for (int dk = -1; dk <= 1; ++dk)
                    {
                        for (int dj = -1; dj <= 1; ++dj)
                        {
                            for (int di = -1; di <= 1; ++di)
                            {
                                if (di == 0 && dj == 0 && dk == 0)
                                {
                                    continue;
                                }

                                if (!twenty_six_conn && std::abs(di) + std::abs(dj) + std::abs(dk) != 1)
                                {
                                    continue;
                                }

                                const int ni = i + di;
                                const int nj = j + dj;
                                const int nk = k + dk;

                                if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
                                {
                                    continue;
                                }

                                const int neighbor_label =
                                    labels[nk][nj][ni];

                                const int neighbor_id =
                                    id(ni, nj, nk);

                                if (neighbor_label < 0)
                                {
                                    touches_anode = true;
                                }
                                else if (neighbor_label > 0)
                                {
                                    touches_cathode = true;
                                }
                                else if (electrolyte[neighbor_id] &&
                                         !visited[neighbor_id])
                                {
                                    visited[neighbor_id] = 1;
                                    q.push({ni, nj, nk});
                                }
                            }
                        }
                    }
                }

                if (touches_anode && touches_cathode)
                {
                    for (const int idx : component)
                    {
                        keep[idx] = 1;
                    }

                    ++kept_components;
                }
                else
                {
                    ++removed_components;
                }
            }
        }
    }

    electrolyte.swap(keep);

    std::cout << "[Full Cell Connectivity] Electrolyte components kept: " << kept_components << ", removed: " << removed_components << "\n";
}



// Half Cell
void Initialize_Geometry::InitializeMesh(const char* meshFile, MPI_Comm comm, int order) {

    myid = mfem::Mpi::WorldRank();

    InitializeGlobalMesh(meshFile);
    InitializeParallelMesh(comm);

    SetupFiniteElementSpace(order);
    SetupParFiniteElementSpace(order);

    AssignGlobalValues();
    MapGlobalToLocal();
    
    particle_labels = GetParticleLabelsFromTiff();

    if (cfg.combine_particle_groups)
    {
        particle_labels.clear();
        particle_labels.push_back(1);

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[Initialize_Geometry] Combining all particle labels into one group.\n";
        }
    }        

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << "[Initialize_Geometry] particle labels found: ";
        for (int lbl : particle_labels) std::cout << lbl << " ";
        std::cout << std::endl;
    }

    HalfCellAMR();

    AllocateHalfCellGeometryFields();
    BuildHalfCellGeometryFields();
    UpdateMeshData();
  
    PrintMeshInfo();

    parallelMesh->SaveAsOne("pmesh");
    MaskFilter->SaveAsOne("MaskFilter.gf");
    MaskFilterPse->SaveAsOne("MaskFilter_pse.gf");
}

void Initialize_Geometry::AllocateHalfCellGeometryFields()
{
    MFEM_VERIFY(parfespace, "Parallel H1 finite element space is not initialized.");

    MaskFilter = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilterPse = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    *MaskFilter = 0.0;
    *MaskFilterPse = 0.0;

    MaskFilters.clear();
    MaskFilters.resize(particle_labels.size());

    for (std::size_t k = 0; k < particle_labels.size(); ++k)
    {
        MaskFilters[k] = std::make_unique<mfem::ParGridFunction>(parfespace.get());
        *MaskFilters[k] = 0.0;
    }
}

void Initialize_Geometry::BuildHalfCellGeometryFields()
{
    MFEM_VERIFY(MaskFilter, "Half-cell total solid mask is not initialized.");
    MFEM_VERIFY(MaskFilterPse, "Half-cell electrolyte mask is not initialized.");

    ComputePDEFilter(*MaskFilter, sim::GeometryPhase::SOLID, CellMode::HALF, cfg.half_electrode);
    ComputePDEFilter(*MaskFilterPse, sim::GeometryPhase::ELECTROLYTE, CellMode::HALF, cfg.half_electrode);

    if (cfg.combine_particle_groups)
    {
        MFEM_VERIFY(MaskFilters.size() == 1, "Expected 1 particle group when combining.");
        *MaskFilters[0] = *MaskFilter;
    }

    else
    {
        const BoundarySide collector_side = (cfg.half_electrode == Electrode::ANODE) ? BoundarySide::WEST : BoundarySide::EAST;

        for (std::size_t k = 0; k < particle_labels.size(); ++k)
        {
            MFEM_VERIFY(MaskFilters[k], "Half-cell particle mask is not allocated.");
            ComputePDEFilterLabel(*MaskFilters[k], particle_labels[k], true, collector_side, CellMode::HALF, cfg.half_electrode);
        }
    }

    if (myid == 0)
    {
        std::cout << "[Initialize_Geometry] Half-cell geometry " << "fields rebuilt on final AMR mesh.\n";
    }
}

void Initialize_Geometry::HalfCellAMR()
{
    if (cfg.amr_levels <= 0)
    {
        return;
    }

    MFEM_VERIFY(parallelMesh, "Parallel mesh is not initialized.");
    MFEM_VERIFY(parfespace, "Parallel H1 space is not initialized.");
    MFEM_VERIFY(parfespace_dg, "Parallel DG space is not initialized.");

    // Temporary fields used only for AMR marking.
    mfem::ParGridFunction temporary_psi(parfespace.get());
    temporary_psi = 0.0;

    ComputePDEFilter(temporary_psi, sim::GeometryPhase::SOLID, CellMode::HALF, cfg.half_electrode);

    const double outer_half_width = 0.495;
    const double size_tolerance = 1.0e-10;

    for (int level = 0; level < cfg.amr_levels; ++level)
    {
        mfem::Array<int> refinement_list;

        const double band_fraction = static_cast<double>(cfg.amr_levels - level) / static_cast<double>(cfg.amr_levels);
        const double half_width = outer_half_width * band_fraction;
        const double psi_lower = 0.5 - half_width;
        const double psi_upper = 0.5 + half_width;

        for (int ei = 0; ei < parallelMesh->GetNE(); ++ei)
        {
            const double element_size = parallelMesh->GetElementSize(ei);

            if (element_size <= cfg.dh * (1.0 + size_tolerance))
            {
                continue;
            }

            mfem::Array<double> values;
            temporary_psi.GetNodalValues(ei, values);

            if (values.Size() == 0)
            {
                continue;
            }

            double element_min = std::numeric_limits<double>::max();
            double element_max = -std::numeric_limits<double>::max();

            for (int j = 0; j < values.Size(); ++j)
            {
                element_min = std::min(element_min, values[j]);
                element_max = std::max(element_max, values[j]);
            }

            const bool intersects_band = element_min < psi_upper && element_max > psi_lower;

            if (intersects_band)
            {
                refinement_list.Append(ei);
            }
        }

        const int local_marked = refinement_list.Size();
        const int local_elements = parallelMesh->GetNE();

        int global_marked = 0;
        int global_elements = 0;

        MPI_Allreduce(&local_marked, &global_marked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_elements, &global_elements, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (myid == 0)
        {
            std::cout << "[AMR] band " << level + 1 << ": psi range = (" << psi_lower << ", " << psi_upper
                << "), marked " << global_marked << " / "
                << global_elements << " elements globally\n";
        }

        if (global_marked == 0)
        {
            if (myid == 0)
            {
                std::cout << "[AMR] No additional elements " << "require refinement.\n";
            }

            break;
        }

        parallelMesh->GeneralRefinement(refinement_list, 1);
        UpdateSpacesAfterAMR();
        temporary_psi.Update();

        ComputePDEFilter(temporary_psi, sim::GeometryPhase::SOLID, CellMode::HALF, cfg.half_electrode);

        PrintAMRMeshInfo(level + 1);
    }
}

void Initialize_Geometry::UpdateSpacesAfterAMR()
{
    MFEM_VERIFY(parfespace, "Parallel H1 space is not initialized.");
    MFEM_VERIFY(parfespace_dg, "Parallel DG space is not initialized.");
    MFEM_VERIFY(pardimfespace_dg, "Parallel vector DG space is not initialized.");

    parfespace->Update();
    parfespace_dg->Update();
    pardimfespace_dg->Update();

    if (Vox)
    {
        Vox->Update();
    }
}

std::vector<std::vector<std::vector<int>>> Initialize_Geometry::MergeMeshes(const char *AnodeMeshFile, const char *CathodeMeshFile)
{
    auto anodeData = ReadTiffFile(AnodeMeshFile);
    auto cathodeData = ReadTiffFile(CathodeMeshFile);

    if (anodeData.empty() ||
        cathodeData.empty() ||
        anodeData[0].empty() ||
        cathodeData[0].empty() ||
        anodeData[0][0].empty() ||
        cathodeData[0][0].empty())
    {
        throw std::runtime_error(
            "Anode or cathode TIFF data is empty.");
    }

    const int anodeNz = static_cast<int>(anodeData.size());
    const int anodeNy = static_cast<int>(anodeData[0].size());
    const int anodeNx = static_cast<int>(anodeData[0][0].size());

    const int cathodeNz = static_cast<int>(cathodeData.size());
    const int cathodeNy = static_cast<int>(cathodeData[0].size());
    const int cathodeNx = static_cast<int>(cathodeData[0][0].size());

    if (anodeNz != cathodeNz)
    {
        throw std::runtime_error(
            "Anode and cathode TIFF files must have "
            "the same number of depth slices.");
    }

    if (anodeNy != cathodeNy)
    {
        throw std::runtime_error(
            "Anode and cathode TIFF files must have "
            "the same number of rows.");
    }

    const int separatorColumns = 0;
    const int mergedNx = anodeNx + separatorColumns + cathodeNx;

    std::vector<std::vector<std::vector<int>>> mergedData( anodeNz,
        std::vector<std::vector<int>>(
            anodeNy,
            std::vector<int>(
                mergedNx,
                0)));

    for (int k = 0; k < anodeNz; ++k)
    {
        for (int j = 0; j < anodeNy; ++j)
        {
            for (int i = 0; i < anodeNx; ++i)
            {
                const int label =
                    anodeData[k][j][i];

                if (label > 0)
                {
                    throw std::runtime_error(
                        "Anode TIFF contains a positive "
                        "particle label. Expected labels <= 0.");
                }

                mergedData[k][j][i] =
                    label;
            }

            const int cathodeStart = anodeNx + separatorColumns;

            for (int i = 0; i < cathodeNx; ++i)
            {
                const int label = cathodeData[k][j][i];

                if (label < 0)
                {
                    throw std::runtime_error(
                        "Cathode TIFF contains a negative "
                        "particle label. Expected labels >= 0.");
                }

                mergedData[k][j][cathodeStart + i] =
                    label;
            }
        }
    }

    if (myid == 0)
    {
        std::cout << "[Full Cell] Signed TIFF files merged.\n";

        std::cout << "  Anode columns:    " << anodeNx << "\n";
        std::cout << "  Separator columns: " << separatorColumns << "\n";
        std::cout << "  Cathode columns:  " << cathodeNx << "\n";
        std::cout << "  Total columns:    " << mergedNx << "\n";
    }

    SaveTiffDataToPGM(mergedData, "full_cell_geometry.pgm");

    return mergedData;
}

// Full Cell
void Initialize_Geometry::InitializeMesh(const char* AnodeMeshFile, const char* CathodeMeshFile, MPI_Comm comm, int order) {

    myid = mfem::Mpi::WorldRank();
    tiffData = MergeMeshes(AnodeMeshFile, CathodeMeshFile);
    
    InitializeGlobalMesh(tiffData);
    InitializeParallelMesh(comm);

    SetupFiniteElementSpace(order);
    SetupParFiniteElementSpace(order);

    AssignGlobalValues();
    MapGlobalToLocal();

    FullCellParticleLabels();

    FullCellAMR();
    AllocateFullCellGeometryFields();
    BuildFullCellGeometryFields();
    UpdateMeshData();

    parallelMesh->SaveAsOne("pmesh");
    MaskFilterAnode->SaveAsOne("MaskFilter_anode.gf");
    MaskFilterCathode->SaveAsOne("MaskFilter_cathode.gf");
    MaskFilterPse->SaveAsOne("MaskFilter_pse.gf");

    if (myid == 0)
    {
        std::cout << "[Initialize_Geometry] " << "Full-cell PDE filters complete.\n";
    }

    PrintMeshInfo();    

}

void Initialize_Geometry::AllocateFullCellGeometryFields()
{
    MFEM_VERIFY(parfespace, "Parallel H1 space is not initialized.");

    MaskFilterAnode = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilterCathode = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilterPse = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    *MaskFilterAnode = 0.0;
    *MaskFilterCathode = 0.0;
    *MaskFilterPse = 0.0;

    MaskFiltersAnode.clear();
    MaskFiltersAnode.resize(anode_particle_labels.size());

    for (std::size_t k = 0; k < anode_particle_labels.size(); ++k)
    {
        MaskFiltersAnode[k] = std::make_unique<mfem::ParGridFunction>(parfespace.get());
        *MaskFiltersAnode[k] = 0.0;
    }

    MaskFiltersCathode.clear();
    MaskFiltersCathode.resize(cathode_particle_labels.size());

    for (std::size_t k = 0; k < cathode_particle_labels.size(); ++k)
    {
        MaskFiltersCathode[k] = std::make_unique<mfem::ParGridFunction>(parfespace.get());
        *MaskFiltersCathode[k] = 0.0;
    }
}

void Initialize_Geometry::BuildFullCellGeometryFields()
{
    ComputePDEFilter(*MaskFilterAnode, sim::GeometryPhase::SOLID, CellMode::FULL, sim::Electrode::ANODE);
    ComputePDEFilter(*MaskFilterCathode, sim::GeometryPhase::SOLID, CellMode::FULL, sim::Electrode::CATHODE);
    ComputePDEFilter(*MaskFilterPse, sim::GeometryPhase::ELECTROLYTE, CellMode::FULL, sim::Electrode::ANODE);

    for (std::size_t k = 0; k < anode_particle_labels.size(); ++k)
    {
        ComputePDEFilterLabel(*MaskFiltersAnode[k], anode_particle_labels[k], false, BoundarySide::WEST, CellMode::FULL, Electrode::ANODE);
    }

    for (std::size_t k = 0; k < cathode_particle_labels.size(); ++k)
    {
        ComputePDEFilterLabel(*MaskFiltersCathode[k], cathode_particle_labels[k], false, BoundarySide::EAST, CellMode::FULL, Electrode::CATHODE);
    }
}

void Initialize_Geometry::FullCellAMR()
{
    if (cfg.amr_levels <= 0)
    {
        return;
    }

    mfem::ParGridFunction temporary_anode(parfespace.get());
    mfem::ParGridFunction temporary_cathode(parfespace.get());
    mfem::ParGridFunction temporary_total(parfespace.get());

    temporary_anode = 0.0;
    temporary_cathode = 0.0;
    temporary_total = 0.0;

    ComputePDEFilter(temporary_anode, sim::GeometryPhase::SOLID, CellMode::FULL, sim::Electrode::ANODE);
    ComputePDEFilter(temporary_cathode, sim::GeometryPhase::SOLID, CellMode::FULL, sim::Electrode::CATHODE);

    temporary_total = temporary_anode;
    temporary_total += temporary_cathode;

    const double outer_half_width = 0.495;
    const double size_tolerance = 1.0e-10;

    for (int level = 0; level < cfg.amr_levels; ++level)
    {
        mfem::Array<int> refinement_list;
        const double band_fraction = static_cast<double>(cfg.amr_levels - level) / static_cast<double>(cfg.amr_levels);
        const double half_width = outer_half_width * band_fraction;
        const double psi_lower = 0.5 - half_width;
        const double psi_upper = 0.5 + half_width;

        for (int ei = 0; ei < parallelMesh->GetNE();++ei)
        {
            const double element_size = parallelMesh->GetElementSize(ei);

            if (element_size <= cfg.dh * (1.0 + size_tolerance))
            {
                continue;
            }

            mfem::Array<double> values;

            temporary_total.GetNodalValues(ei, values);

            if (values.Size() == 0)
            {
                continue;
            }

            double element_min = std::numeric_limits<double>::max();
            double element_max = -std::numeric_limits<double>::max();

            for (int j = 0; j < values.Size(); ++j)
            {
                element_min = std::min(element_min, values[j]);
                element_max = std::max(element_max, values[j]);
            }

            if (element_min < psi_upper && element_max > psi_lower)
            {
                refinement_list.Append(ei);
            }
        }

        const int local_marked = refinement_list.Size();
        int global_marked = 0;

        MPI_Allreduce(&local_marked, &global_marked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (global_marked == 0)
        {
            break;
        }

        parallelMesh->GeneralRefinement(refinement_list, 1);

        UpdateSpacesAfterAMR();

        temporary_anode.Update();
        temporary_cathode.Update();
        temporary_total.Update();

        ComputePDEFilter(temporary_anode, sim::GeometryPhase::SOLID, sim::CellMode::FULL, sim::Electrode::ANODE);
        ComputePDEFilter(temporary_cathode, sim::GeometryPhase::SOLID, sim::CellMode::FULL, sim::Electrode::CATHODE);

        temporary_total = temporary_anode;
        temporary_total += temporary_cathode;

        PrintAMRMeshInfo(level + 1);
    }
}



std::vector<int> Initialize_Geometry::GetParticleLabelsFromTiff() const
{
    std::set<int> labels;

    for (const auto &slice : tiffData)
    for (const auto &row   : slice)
    for (int v : row)
    {
        if (v != 0) { labels.insert(v); } // 0 is electrolyte
    }

    return std::vector<int>(labels.begin(), labels.end());
}

void Initialize_Geometry::FullCellParticleLabels()
{
    std::set<int> anodeLabels;
    std::set<int> cathodeLabels;

    for (const auto &slice : tiffData)
    {
        for (const auto &row : slice)
        {
            for (const int label : row)
            {
                if (label < 0)
                {
                    anodeLabels.insert(label);
                }
                else if (label > 0)
                {
                    cathodeLabels.insert(label);
                }
            }
        }
    }

    anode_particle_labels.assign(anodeLabels.begin(), anodeLabels.end());
    cathode_particle_labels.assign(cathodeLabels.begin(), cathodeLabels.end());

    std::sort(anode_particle_labels.begin(), anode_particle_labels.end(), [](const int lhs, const int rhs)
        {
            return std::abs(lhs) < std::abs(rhs);
        }
    );

    if (cfg.combine_particle_groups)
    {
        if (!anode_particle_labels.empty())
        {
            anode_particle_labels.clear();
            anode_particle_labels.push_back(-1);
        }

        if (!cathode_particle_labels.empty())
        {
            cathode_particle_labels.clear();
            cathode_particle_labels.push_back(1);
        }

        if (myid == 0)
        {
            std::cout << "[Initialize_Geometry] " << "Combining particles separately within " << "each full-cell electrode.\n";
        }
    }

    if (myid == 0)
    {
        std::cout << "[Initialize_Geometry] " << "Anode particle labels: ";

        for (const int label : anode_particle_labels)
        {
            std::cout << label << " ";
        }

        std::cout << "\n[Initialize_Geometry] " << "Cathode particle labels: ";

        for (const int label : cathode_particle_labels)
        {
            std::cout << label << " ";
        }

        std::cout << "\n";
    }
}

// Function to initialize the global mesh using a .tif or .mesh file
void Initialize_Geometry::InitializeGlobalMesh(const char* meshFile) {
    std::string meshFileStr(meshFile);  // Convert to std::string
    std::string fileExtension = meshFileStr.substr(meshFileStr.find_last_of(".") + 1);

    if (fileExtension == "tif") {
        if (mfem::Mpi::WorldRank() == 0) { std::cout << "Creating global mesh using .tif file" << std::endl; }
        tiffData = ReadTiffFile(meshFile); // read voxel data from tiff file
        globalMesh = CreateGlobalMeshFromTiffData(tiffData); // generate mesh from voxel data
    } 
    else {
        throw std::invalid_argument("Unsupported file format. Only .tif is allowed.");
    }

    // ensure mesh supports non-conforming elements for adaptive refinement
    globalMesh->EnsureNCMesh(true);

    int e = 0;
    mfem::Array<int> vert_ids;
    globalMesh->GetElementVertices(e, vert_ids);

    mfem::Vector v0(globalMesh->GetVertex(vert_ids[0]), globalMesh->SpaceDimension());
    mfem::Vector v1(globalMesh->GetVertex(vert_ids[1]), globalMesh->SpaceDimension());

    double dh1 = v0.DistanceTo(v1);
    if (mfem::Mpi::WorldRank() == 0) { std::cout << "Distributed element size dh = " << dh1 << std::endl;}

    MPI_Barrier(MPI_COMM_WORLD);

}

void Initialize_Geometry::InitializeGlobalMesh( const std::vector<std::vector<std::vector<int>>> &voxelData)
{
    if (voxelData.empty() || voxelData[0].empty() || voxelData[0][0].empty())
    {
        throw std::invalid_argument("InitializeGlobalMesh: voxel data is empty.");
    }

    if (myid == 0)
    {
        std::cout << "Creating global mesh from voxel data " << "already stored in memory.\n";
    }

    // Copy merged data into the class member.
    tiffData = voxelData;

    globalMesh = CreateGlobalMeshFromTiffData(tiffData);
    globalMesh->EnsureNCMesh(true);

    if (globalMesh->GetNE() <= 0)
    {
        throw std::runtime_error("InitializeGlobalMesh: generated mesh has no elements.");
    }

    mfem::Array<int> vertexIds;
    globalMesh->GetElementVertices(0, vertexIds);

    mfem::Vector vertex0(globalMesh->GetVertex(vertexIds[0]), globalMesh->SpaceDimension());
    mfem::Vector vertex1(globalMesh->GetVertex(vertexIds[1]), globalMesh->SpaceDimension());

    const double elementSize = vertex0.DistanceTo(vertex1);

    if (myid == 0)
    {
        std::cout << "Element size dh = " << elementSize << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// Function to initialize the parallel mesh
void Initialize_Geometry::InitializeParallelMesh(MPI_Comm comm) {
    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before creating a parallel mesh.");
    }
    parallelMesh = std::make_shared<mfem::ParMesh>(comm, *globalMesh);
    parallelMesh->SaveAsOne("pmesh");

}

// Function to set up the finite element space on global mesh
void Initialize_Geometry::SetupFiniteElementSpace(int order) {
    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before setting up FE space.");
    }

    gfec = std::make_unique<mfem::H1_FECollection>(order, globalMesh->Dimension());
    globalfespace = std::make_shared<mfem::FiniteElementSpace>(globalMesh.get(), gfec.get());
}

// Function to set up finite element space on parallel mesh
void Initialize_Geometry::SetupParFiniteElementSpace(int order) {
    if (!parallelMesh) {
        throw std::runtime_error("Parallel mesh must be initialized before setting up FE space.");
    }
    
    pfec = std::make_unique<mfem::H1_FECollection>(order, parallelMesh->Dimension());
    parfespace = std::make_shared<mfem::ParFiniteElementSpace>(parallelMesh.get(), pfec.get());

    this->pfec_dg = std::make_unique<mfem::DG_FECollection>(order, this->parallelMesh->Dimension(), mfem::BasisType::GaussLobatto);
    this->parfespace_dg = std::make_shared<mfem::ParFiniteElementSpace>(this->parallelMesh.get(), this->pfec_dg.get());
    this->pardimfespace_dg = std::make_shared<mfem::ParFiniteElementSpace>(this->parallelMesh.get(), this->pfec_dg.get(), this->parallelMesh->Dimension(), mfem::Ordering::byNODES);
}

void Initialize_Geometry::AssignGlobalValues()
{
    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << "Reading TIFF file for voxel data\n";
    }

    if (!globalfespace)
    {
        throw std::runtime_error("Global finite element space must be initialized before assigning global values.");
    }

    gVox = std::make_unique<mfem::GridFunction>(globalfespace.get());

    const int nz = static_cast<int>(tiffData.size());
    const int ny = static_cast<int>(tiffData[0].size());
    const int nx = static_cast<int>(tiffData[0][0].size());

    const int coarsen = cfg.coarsen_factor;

    const int ex = (nx - 1) / coarsen;
    const int ey = (ny - 1) / coarsen;

    const int vx = ex + 1;
    const int vy = ey + 1;

    *gVox = 0.0;

    if (nz == 1)
    {
        for (int j = 0; j < vy; ++j)
        {
            for (int i = 0; i < vx; ++i)
            {
                const int ii = coarsen * i;
                const int jj = coarsen * j;

                const int idx = i + vx * j;

                (*gVox)[idx] = tiffData[0][jj][ii];
            }
        }
    }
    else
    {
        const int ez = (nz - 1) / coarsen;
        const int vz = ez + 1;

        for (int k = 0; k < vz; ++k)
        {
            for (int j = 0; j < vy; ++j)
            {
                for (int i = 0; i < vx; ++i)
                {
                    const int ii = coarsen * i;
                    const int jj = coarsen * j;
                    const int kk = coarsen * k;

                    const int idx =
                        i + vx * (j + vy * k);

                    (*gVox)[idx] = tiffData[kk][jj][ii];
                }
            }
        }
    }
}

void Initialize_Geometry::MapGlobalToLocal() {
    
    if (!parallelMesh) {
        throw std::runtime_error("Parallel mesh must be initialized before calculating element volumes.");
    }

    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before setting up FE space.");
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
    nV = parallelMesh->GetNV();        // number of vertices
    nE = parallelMesh->GetNE();        // number of elements
    nC = pow(2, parallelMesh->Dimension());  // number of corner vertices

    parallelMesh->GetGlobalElementIndices(E_L2G);

    gVTX.SetSize(nC);
    VTX.SetSize(nC);

    if (mfem::Mpi::WorldRank() == 0) // only print on rank 0
    {cout << "Reading .tif file for mapping global to local grid function" << endl;}

    Vox = std::make_unique<mfem::ParGridFunction>(parfespace.get()); // used in Vox code

    // Iterate over elements and map global to local
    for (ei = 0; ei < nE; ei++) {
        gei = E_L2G[ei];

        globalMesh->GetElementVertices(gei, gVTX);
        parallelMesh->GetElementVertices(ei, VTX);

        for (int vi = 0; vi < nC; vi++) {                            // used in Vox code
            (*this->Vox)(VTX[vi]) = (*this->gVox)(gVTX[vi]);         // used in Vox code
        }   
        
    }
}

// Reading .tif file and returning voxel data
std::vector<std::vector<std::vector<int>>> Initialize_Geometry::ReadTiffFile(const char* meshFile) {

	if (mfem::Mpi::WorldRank() == 0) { std::cout << "reading tiff file" << std::endl; }
	Constraints args;
    
    args.Depth_begin = cfg.depth_begin;
    args.Depth_end   = cfg.depth_end;

    args.Row_begin = cfg.row_begin;
    args.Row_end   = cfg.row_end;

    args.Column_begin = cfg.column_begin;
    args.Column_end   = cfg.column_end;

	TIFFReader reader(meshFile,args);
	reader.readinfo();
	std::vector<std::vector<std::vector<int>>> tiffData;
	tiffData = reader.getImageData();

    return tiffData;
}

std::unique_ptr<mfem::Mesh> Initialize_Geometry::CreateGlobalMeshFromTiffData(const std::vector<std::vector<std::vector<int>>>& tiffData)
{
    const int nz = static_cast<int>(tiffData.size());
    const int ny = static_cast<int>(tiffData[0].size());
    const int nx = static_cast<int>(tiffData[0][0].size());

    const int coarsen = cfg.coarsen_factor;

    const int ex = (nx - 1) / coarsen;
    const int ey = (ny - 1) / coarsen;
    const int ez = (nz == 1) ? 1 : (nz - 1) / coarsen;

    // Preserve the physical size of the selected TIFF region.
    const double sx = (nx - 1) * cfg.dh;
    const double sy = (ny - 1) * cfg.dh;
    const double sz = (nz == 1) ? cfg.dh : (nz - 1) * cfg.dh;

    const bool generate_edges = false;
    const bool sfc_ordering = false;

    std::unique_ptr<mfem::Mesh> mesh;

    if (nz == 1)
    {
        mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(ex, ey, mfem::Element::QUADRILATERAL, generate_edges, sx, sy, sfc_ordering));
    }
    else
    {
        mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(ex, ey, ez, mfem::Element::HEXAHEDRON, sx, sy, sz, sfc_ordering));
    }

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << "Initial mesh dimensions:\n";
        std::cout << "  elements x = " << ex << "\n";
        std::cout << "  elements y = " << ey << "\n";

        if (nz > 1)
        {
            std::cout << "  elements z = `" << ez << "\n";
        }

        std::cout << "  element size after coarsening = " << cfg.dh * coarsen << "\n";
    }

    return mesh;
}

void Initialize_Geometry::PrintMeshInfo() {
    
    if (!parallelMesh) {
        std::cout << "Parallel mesh not initialized.\n";
        return;
    }

}

void Initialize_Geometry::SaveTiffDataToPGM(const std::vector<std::vector<std::vector<int>>> &data, const std::string &filename)
{
    if (data.empty() || data[0].empty() || data[0][0].empty()) {
        std::cerr << "SaveTiffDataToPGM: empty data\n";
        return;
    }

    const auto &img = data[0];
    const int height = static_cast<int>(img.size());
    const int width  = static_cast<int>(img[0].size());

    int minLabel = std::numeric_limits<int>::max();
    int maxLabel = std::numeric_limits<int>::lowest();

    for (const auto &row : img)
    {
        for (const int label : row)
        {
            minLabel = std::min(minLabel, label);
            maxLabel = std::max(maxLabel, label);
        }
    }

    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Could not open file for writing: " << filename << "\n";
        return;
    }

    out << "P5\n" << width << " " << height << "\n255\n";

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            const int label = img[j][i];

            unsigned char val = 0;
            if (maxLabel > minLabel) {
                val = static_cast<unsigned char>(std::round(255.0 * (label - minLabel) / (maxLabel - minLabel)));
            }

            out.write(reinterpret_cast<char*>(&val), 1);
        }
    }

    out.close();

    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "Saved PGM to " << filename
                  << " using label range " << minLabel << "-" << maxLabel << "\n";
    }
}

void Initialize_Geometry::ApplyPDEFilterToMask(const std::vector<uint8_t>& mask, int nx, int ny, int nz, mfem::ParGridFunction& filt_gf)
{
    mfem::ParGridFunction ls_coeff_dg(parfespace_dg.get());
    mfem::ParGridFunction filt_dg(parfespace_dg.get());

    ls_coeff_dg = 0.0;
    filt_dg = 0.0;

    struct MaskCoefficient : public mfem::Coefficient
    {
        int nx;
        int ny;
        int nz;
        int dim;

        double x0;
        double y0;
        double z0;

        double dx;
        double dy;
        double dz;

        const std::vector<uint8_t>* mask;

        MaskCoefficient(int nx_, int ny_, int nz_, mfem::ParMesh& mesh, const std::vector<uint8_t>& mask_, double fallback_spacing)
            : nx(nx_), ny(ny_), nz(nz_), dim(mesh.Dimension()), mask(&mask_)
        {
            mfem::Vector bb_min;
            mfem::Vector bb_max;

            mesh.GetBoundingBox(bb_min, bb_max);

            x0 = bb_min(0);
            y0 = bb_min(1);
            z0 = (dim == 3) ? bb_min(2) : 0.0;

            dx = (nx > 1) ? (bb_max(0) - bb_min(0)) / (nx - 1) : fallback_spacing;
            dy = (ny > 1) ? (bb_max(1) - bb_min(1)) / (ny - 1) : fallback_spacing;
            dz = (dim == 3 && nz > 1) ? (bb_max(2) - bb_min(2)) / (nz - 1) : fallback_spacing;
        }

        double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override
        {
            mfem::Vector x;
            T.Transform(ip, x);

            int i = static_cast<int>(std::floor((x(0) - x0) / dx + 0.5));
            int j = static_cast<int>(std::floor((x(1) - y0) / dy + 0.5));
            int k = 0;

            if (dim == 3)
            {
                k = static_cast<int>(std::floor((x(2) - z0) / dz + 0.5));
            }

            i = std::max(0, std::min(nx - 1, i));
            j = std::max(0, std::min(ny - 1, j));
            k = std::max(0, std::min(nz - 1, k));

            const int idx = i + nx * j + nx * ny * k;

            return (*mask)[idx] ? 1.0 : -1.0;
        }
    };

    MaskCoefficient mask_coefficient(nx, ny, nz, *parallelMesh, mask, cfg.dh);

    ls_coeff_dg.ProjectCoefficient(mask_coefficient);
    const double filter_weight = 3.0 * cfg.dh;

    mfem::common::PDEFilter filter(*parallelMesh, filter_weight);
    filter.Filter(ls_coeff_dg, filt_dg);

    for (int i = 0; i < filt_dg.Size(); ++i)
    {
        filt_dg(i) = 0.5 * (filt_dg(i) + 1.0);
    }

    mfem::GridFunctionCoefficient filtered_coefficient(&filt_dg);

    filt_gf.ProjectCoefficient(filtered_coefficient);

    mfem::Vector true_values;
    filt_gf.GetTrueDofs(true_values);
    filt_gf.SetFromTrueDofs(true_values);
}


void Initialize_Geometry::ComputePDEFilter(mfem::ParGridFunction &filt_gf, sim::GeometryPhase phase, sim::CellMode cell_mode, sim::Electrode electrode)

{
    MFEM_VERIFY(parallelMesh, "parallelMesh is not initialized.");
    MFEM_VERIFY(parfespace, "parfespace is not initialized.");
    MFEM_VERIFY(filt_gf.ParFESpace() == parfespace.get(), "filt_gf must be on parfespace.");
    MFEM_VERIFY(parfespace_dg, "parfespace_dg is not initialized.");
    MFEM_VERIFY(parallelMesh->Dimension() == 2 || parallelMesh->Dimension() == 3, "ComputePDEFilter: mesh must be 2D or 3D.");

    // TIFF sizes
    const int nz = (int)tiffData.size();
    const int ny = (int)tiffData[0].size();
    const int nx = (int)tiffData[0][0].size();

    const int rank = mfem::Mpi::WorldRank();
    const bool eight_conn = false;        // 2D: 4/8
    const bool twenty_six = false;        // 3D: 6/26  (set true if desired)

    std::vector<uint8_t> fg(nx*ny*nz, 0);

    // Boundary Rules: [0] = west, [1] = east, [2] = south, [3] = north, [4] = bottom, [5] = top
    if (rank == 0)
    {
        for (int k = 0; k < nz; ++k)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int i = 0; i < nx; ++i)
                {
                    const int idx = i + nx*j + nx*ny*k;
                    const int label = tiffData[k][j][i];

                    if (phase == GeometryPhase::SOLID)
                    {
                        // Electrode mask.
                        if (cell_mode == CellMode::HALF)
                        {
                            fg[idx] = (label > 0) ? 1 : 0;
                        }
                        else if (electrode == Electrode::ANODE)
                        {
                            fg[idx] = (label < 0) ? 1 : 0;
                        }
                        else if (electrode == Electrode::CATHODE)
                        {
                            fg[idx] = (label > 0) ? 1 : 0;
                        }
                        else
                        {
                            MFEM_ABORT("ComputePDEFilter: full-cell solid filter requires ANODE or CATHODE.");
                        }
                    }
                    else if (phase == GeometryPhase::ELECTROLYTE)
                    {
                        fg[idx] = (label == 0) ? 1 : 0;
                    }
                    else
                    {
                        MFEM_ABORT("ComputePDEFilter: phase must be SOLID or ELECTROLYTE.");
                    }
                }
            }
        }

        if (nz == 1)
        {
            if (phase == GeometryPhase::SOLID)
            {
                const BoundarySide collectorSide = (electrode == Electrode::ANODE) ? BoundarySide::WEST : BoundarySide::EAST;
                KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, collectorSide);
            }
            else
            {
                if (cell_mode == CellMode::FULL)
                {
                    KeepOnlyElectrolyteTouchingBothElectrodes_2D(fg, tiffData, nx, ny, eight_conn);
                }
                else
                {
                    const BoundarySide electrolyteSide = (electrode == Electrode::ANODE) ? BoundarySide::EAST : BoundarySide::WEST;
                    KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, electrolyteSide);

                }
            }
        }
        else
        {
            if (phase == GeometryPhase::SOLID)
            {
                const BoundarySide collectorFace = (electrode == Electrode::ANODE) ? BoundarySide::WEST : BoundarySide::EAST;
                KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, false, collectorFace);
            }
            else
            {
                if (cell_mode == CellMode::FULL)
                {
                    KeepOnlyElectrolyteTouchingBothElectrodes_3D(fg, tiffData, nx, ny, nz, twenty_six);
                }
                else
                {
                    const BoundarySide electrolyteFace = (electrode == Electrode::ANODE) ? BoundarySide::EAST : BoundarySide::WEST;
                    KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, false, electrolyteFace);
                }
            }
        }
    }

    MPI_Bcast(fg.data(), (int)fg.size(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    ApplyPDEFilterToMask(fg, nx, ny, nz, filt_gf);
}

void Initialize_Geometry::ComputePDEFilterLabel(mfem::ParGridFunction &filt_gf, int target_label,
                                                bool keep_boundary_connected, BoundarySide seed_side_or_face, sim::CellMode cell_mode, sim::Electrode electrode)
{
    MFEM_VERIFY(parallelMesh, "parallelMesh is not initialized.");
    MFEM_VERIFY(parfespace, "parfespace is not initialized.");
    MFEM_VERIFY(filt_gf.ParFESpace() == parfespace.get(), "filt_gf must be on parfespace.");
    MFEM_VERIFY(parfespace_dg, "parfespace_dg is not initialized.");
    MFEM_VERIFY(parallelMesh->Dimension() == 2 || parallelMesh->Dimension() == 3, "ComputePDEFilterLabel: mesh must be 2D or 3D.");

    const int nz = (int)tiffData.size();
    const int ny = (int)tiffData[0].size();
    const int nx = (int)tiffData[0][0].size();

    const int rank = mfem::Mpi::WorldRank();
    const bool eight_conn = false;
    const bool twenty_six = false;

    std::vector<uint8_t> fg(nx * ny * nz, 0);

    if (rank == 0)
    {
        // Mask containing every solid particle, regardless of label.
        std::vector<uint8_t> connected_solid(nx * ny * nz, 0);

        for (int k = 0; k < nz; ++k)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int i = 0; i < nx; ++i)
                {
                    const int idx = i + nx*j + nx*ny*k;
                    const int label = tiffData[k][j][i];

                    // Half-cell: any nonzero/positive particle label is solid.
                    if (cell_mode == CellMode::HALF)
                    {
                        connected_solid[idx] = (label > 0) ? 1 : 0;
                    }
                    else if (electrode == Electrode::ANODE)
                    {
                        connected_solid[idx] = (label < 0) ? 1 : 0;
                    }
                    else if (electrode == Electrode::CATHODE)
                    {
                        connected_solid[idx] = (label > 0) ? 1 : 0;
                    }
                }
            }
        }

        // Find the entire particle network connected to the collector.
        if (keep_boundary_connected)
        {
            if (nz == 1)
            {
                KeepOnlyConnectedToBoundary_2D(connected_solid, nx, ny, eight_conn, false, seed_side_or_face);
            }
            else
            {
                KeepOnlyConnectedToBoundary_3D(connected_solid, nx, ny, nz, twenty_six, false, seed_side_or_face);
            }
        }

        for (int k = 0; k < nz; ++k)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int i = 0; i < nx; ++i)
                {
                    const int idx = i + nx*j + nx*ny*k;
                    const int label = tiffData[k][j][i];

                    bool belongs_to_group = false;

                    if (!cfg.combine_particle_groups)
                    {
                        belongs_to_group = (label == target_label);
                    }
                    else if (cell_mode == CellMode::HALF)
                    {
                        belongs_to_group = (label > 0);
                    }
                    else if (electrode == Electrode::ANODE)
                    {
                        belongs_to_group = (label < 0);
                    }
                    else if (electrode == Electrode::CATHODE)
                    {
                        belongs_to_group = (label > 0);
                    }

                    fg[idx] =
                        (belongs_to_group && connected_solid[idx]) ? 1 : 0;
                }
            }
        }
    }

    MPI_Bcast(fg.data(), (int)fg.size(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    ApplyPDEFilterToMask(fg, nx, ny, nz, filt_gf);
}

void Initialize_Geometry::UpdateMeshData()
{
    MFEM_VERIFY(parallelMesh, "Parallel mesh is not initialized.");

    nV = parallelMesh->GetNV();
    nE = parallelMesh->GetNE();

    if (nE > 0)
    {
        nC = parallelMesh->GetElement(0)->GetNVertices();
    }
    else
    {
        nC = 0;
    }
}

void Initialize_Geometry::PrintAMRMeshInfo(int level) const
{
    double local_hmin = std::numeric_limits<double>::max();
    double local_hmax = 0.0;

    for (int ei = 0; ei < parallelMesh->GetNE(); ++ei)
    {
        const double h = parallelMesh->GetElementSize(ei);
        local_hmin = std::min(local_hmin, h);
        local_hmax = std::max(local_hmax, h);
    }

    double global_hmin = 0.0;
    double global_hmax = 0.0;

    MPI_Allreduce(&local_hmin, &global_hmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_hmax, &global_hmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const int local_elements = parallelMesh->GetNE();
    int global_elements = 0;

    MPI_Allreduce(&local_elements, &global_elements, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (myid == 0)
    {
        std::cout << "[AMR] level " << level
            << " complete: h_min = " << global_hmin << ", h_max = " << global_hmax
            << ", total elements = " << global_elements << "\n";
    }
}
