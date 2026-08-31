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


// Constructor
Initialize_Geometry::Initialize_Geometry(const SimulationConfig& cfg)
    : cfg(cfg)
{}

// Destructor
Initialize_Geometry::~Initialize_Geometry() {}


static void KeepOnlyConnectedToBoundary_2D(std::vector<uint8_t> &solid, int nx, int ny, bool eight_conn,
                                          bool seed_all_boundaries = true, int seed_side = -1)
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
    if (seed_all_boundaries || seed_side == -1)
    {
        for (int i=0;i<nx;i++){ push(i,0); push(i,ny-1); }
        for (int j=0;j<ny;j++){ push(0,j); push(nx-1,j); }
    }
    else
    {
        if (seed_side == 0) for (int j=0;j<ny;j++) push(0,j);         // left
        if (seed_side == 1) for (int j=0;j<ny;j++) push(nx-1,j);      // right
        if (seed_side == 2) for (int i=0;i<nx;i++) push(i,0);         // bottom
        if (seed_side == 3) for (int i=0;i<nx;i++) push(i,ny-1);      // top
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

static void KeepOnlyConnectedToBoundary_3D(std::vector<uint8_t> &solid,
                                          int nx, int ny, int nz,
                                          bool twenty_six_conn,
                                          bool seed_all_boundaries = true,
                                          int seed_face = -1)
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

    auto seed_face_fn = [&](int face){
        if (face==0) for (int k=0;k<nz;k++) for (int j=0;j<ny;j++) push(0,j,k);         // xmin
        if (face==1) for (int k=0;k<nz;k++) for (int j=0;j<ny;j++) push(nx-1,j,k);      // xmax
        if (face==2) for (int k=0;k<nz;k++) for (int i=0;i<nx;i++) push(i,0,k);         // ymin
        if (face==3) for (int k=0;k<nz;k++) for (int i=0;i<nx;i++) push(i,ny-1,k);      // ymax
        if (face==4) for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) push(i,j,0);         // zmin
        if (face==5) for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) push(i,j,nz-1);      // zmax
    };

    if (seed_all_boundaries || seed_face == -1)
    {
        for (int face=0; face<6; face++) seed_face_fn(face);
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

static void KeepOnlyElectrolyteTouchingBothElectrodes_2D(
    std::vector<uint8_t> &electrolyte,
    const std::vector<std::vector<std::vector<int>>> &labels,
    int nx,
    int ny,
    bool eight_conn)
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

    std::cout
        << "[Full Cell Connectivity] Electrolyte components kept: "
        << kept_components
        << ", removed: "
        << removed_components
        << "\n";
}

static void KeepOnlyElectrolyteTouchingBothElectrodes_3D(
    std::vector<uint8_t> &electrolyte,
    const std::vector<std::vector<std::vector<int>>> &labels,
    int nx,
    int ny,
    int nz,
    bool twenty_six_conn)
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

                                if (!twenty_six_conn &&
                                    std::abs(di) +
                                    std::abs(dj) +
                                    std::abs(dk) != 1)
                                {
                                    continue;
                                }

                                const int ni = i + di;
                                const int nj = j + dj;
                                const int nk = k + dk;

                                if (ni < 0 || ni >= nx ||
                                    nj < 0 || nj >= ny ||
                                    nk < 0 || nk >= nz)
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

    std::cout
        << "[Full Cell Connectivity] Electrolyte components kept: "
        << kept_components
        << ", removed: "
        << removed_components
        << "\n";
}



// Half Cell
void Initialize_Geometry::InitializeMesh(const char* meshFile, MPI_Comm comm, int order) {

    myid = mfem::Mpi::WorldRank();

    // Initialize the global mesh
    InitializeGlobalMesh(meshFile);

    // Initialize the parallel mesh
    InitializeParallelMesh(comm);

    // Set up the finite element space
    SetupFiniteElementSpace(order);

    // Set up the parallel finite element space
    SetupParFiniteElementSpace(order);

    // Assign the global values
    AssignGlobalValues();

    // Map the global values to the local
    MapGlobalToLocal();
    
    // std::string meshFileStr(meshFile);


    // if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "tif")
    // {
        // distMask       = std::make_unique<mfem::ParGridFunction>(parfespace.get());
        // distMaskSigned = std::make_unique<mfem::ParGridFunction>(parfespace.get());

        // MaskFilter    = std::make_unique<mfem::ParGridFunction>(parfespace.get());   // total solid
        // MaskFilterPse = std::make_unique<mfem::ParGridFunction>(parfespace.get());   // electrolyte

        // // Keep your old total-solid and electrolyte filters
        // ComputePDEFilter(*distMask, *MaskFilter,    /*mode=*/0, sim::CellMode::HALF, cfg.half_electrode);
        // ComputePDEFilter(*distMask, *MaskFilterPse, /*mode=*/1, sim::CellMode::HALF, cfg.half_electrode);

        // discover particle labels automatically from TIFF
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

        // // allocate one filtered mask per particle label
        // MaskFilters.clear();
        // MaskFilters.resize(particle_labels.size());

        // for (int k = 0; k < (int)particle_labels.size(); ++k)
        // {
        //     MaskFilters[k] = std::make_unique<mfem::ParGridFunction>(parfespace.get());
        //     ComputePDEFilterLabel(*distMask, *MaskFilters[k], particle_labels[k], false, -1, CellMode::HALF, cfg.half_electrode);

        //     std::ostringstream name;
        //     name << "MaskFilter_label_" << particle_labels[k] << ".gf";
        //     // MaskFilters[k]->SaveAsOne(name.str().c_str());
        // }

        // if (mfem::Mpi::WorldRank() == 0) {
        //     std::cout << "ComputePDEFilter done.\n";
        // }
        
        PrintMeshInfo();

        parallelMesh->SaveAsOne("pmesh");

        MaskFilter->SaveAsOne("MaskFilter.gf");
        MaskFilterPse->SaveAsOne("MaskFilter_pse.gf");

        // if (mfem::Mpi::WorldRank() == 0) {
        //     std::cout << "ComputePDEFilter done.\n";
        // }

    // }

    // Print out information relative to the mesh
    // PrintMeshInfo();

    // globalMesh->Save("gmesh");


}

void Initialize_Geometry::AllocateHalfCellGeometryFields()
{
    MFEM_VERIFY(parfespace, "Parallel H1 finite element space is not initialized.");

    distMask = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    distMaskSigned = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilter = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilterPse = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    *distMask = 0.0;
    *distMaskSigned = 0.0;
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
    MFEM_VERIFY(distMask, "Half-cell filter workspace is not initialized.");
    MFEM_VERIFY(MaskFilter, "Half-cell total solid mask is not initialized.");
    MFEM_VERIFY(MaskFilterPse, "Half-cell electrolyte mask is not initialized.");

    ComputePDEFilter(*distMask, *MaskFilter, 0, CellMode::HALF, cfg.half_electrode);
    ComputePDEFilter(*distMask, *MaskFilterPse, 1, CellMode::HALF, cfg.half_electrode);

    if (cfg.combine_particle_groups)
    {
        MFEM_VERIFY(MaskFilters.size() == 1, "Expected 1 particle group when combining.");
        *MaskFilters[0] = *MaskFilter;
    }
    // if (cfg.combine_particle_groups)
    // {
    //     MFEM_VERIFY(
    //         MaskFilters.size() == 1,
    //         "Expected 1 particle group when combining.");

    //     ComputePDEFilterLabel(
    //         *distMask,
    //         *MaskFilters[0],
    //         1,
    //         false,
    //         -1,
    //         CellMode::HALF,
    //         cfg.half_electrode);
    // }
    else
    {
        const int collector_side = (cfg.half_electrode == Electrode::ANODE) ? 0 : 1;

        for (std::size_t k = 0; k < particle_labels.size(); ++k)
        {
            MFEM_VERIFY(MaskFilters[k], "Half-cell particle mask is not allocated.");
            // const int collector_side = (cfg.half_electrode == Electrode::ANODE) ? 0 : 1;

            ComputePDEFilterLabel(*distMask, *MaskFilters[k], particle_labels[k], true, collector_side, CellMode::HALF, cfg.half_electrode);
            // ComputePDEFilterLabel(*distMask, *MaskFilters[k], particle_labels[k], false, -1, CellMode::HALF, cfg.half_electrode);
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
    mfem::ParGridFunction temporary_dist(parfespace.get());
    mfem::ParGridFunction temporary_psi(parfespace.get());

    temporary_dist = 0.0;
    temporary_psi = 0.0;

    ComputePDEFilter(temporary_dist, temporary_psi, 0, CellMode::HALF, cfg.half_electrode);

    // {
    // std::ostringstream name;
    // name << "psi_amr_level_0.gf";
    // temporary_psi.SaveAsOne(name.str().c_str());
    // }

    // {
    // std::ostringstream meshname;
    // meshname << "pmesh_before_amr";
    // parallelMesh->SaveAsOne(meshname.str().c_str());
    // }

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

        temporary_dist.Update();
        temporary_psi.Update();

        ComputePDEFilter(temporary_dist, temporary_psi, 0, CellMode::HALF, cfg.half_electrode);

        // {
        //     std::ostringstream name;
        //     name << "psi_amr_level_" << (level + 1) << ".gf";
        //     temporary_psi.SaveAsOne(name.str().c_str());
        // }

        // {
        //     std::ostringstream meshname;
        //     meshname << "pmesh_amr_level_" << (level + 1);
        //     parallelMesh->SaveAsOne(meshname.str().c_str());
        // }
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

    /*
     * This initially joins the two electrode geometries
     * directly. Set this to a positive integer later if
     * you want an explicit electrolyte-only separator.
     */
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
            // Copy the anode to the left side.
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

            /*
             * Separator entries remain zero because
             * mergedData was initialized to zero.
             */

            const int cathodeStart = anodeNx + separatorColumns;

            // Copy the cathode to the right side.
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
    
    // Initialize the global mesh
    InitializeGlobalMesh(tiffData);

    // Initialize the parallel mesh
    InitializeParallelMesh(comm);

    // Set up the finite element space
    SetupFiniteElementSpace(order);

    // Set up the parallel finite element space
    SetupParFiniteElementSpace(order);

    // Assign the global values
    AssignGlobalValues();

    // Map the global values to the local
    MapGlobalToLocal();

    // // General filter workspace.
    // distMask = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    // distMaskSigned = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    // // Full-cell fields.
    // MaskFilterAnode = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    // MaskFilterCathode = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    // MaskFilterPse = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    // ComputePDEFilter(*distMask, *MaskFilterAnode, 0, CellMode::FULL, Electrode::ANODE);
    // ComputePDEFilter(*distMask, *MaskFilterCathode, 0, CellMode::FULL, Electrode::CATHODE);

    // // Zero labels: electrolyte.
    // // Electrode argument is ignored when phase_mode == 1.
    // ComputePDEFilter(*distMask, *MaskFilterPse, 1, CellMode::FULL, Electrode::ANODE);

    // Discover negative and positive labels separately.
    FullCellParticleLabels();

    FullCellAMR();
    AllocateFullCellGeometryFields();
    BuildFullCellGeometryFields();
    UpdateMeshData();

    // // -------------------------------------------------
    // // Anode particle masks
    // // -------------------------------------------------

    // MaskFiltersAnode.clear();
    // MaskFiltersAnode.resize(anode_particle_labels.size());

    // for (int p = 0;
    //     p < static_cast<int>(anode_particle_labels.size());
    //     ++p)
    // {
    //     MaskFiltersAnode[p] = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    //     ComputePDEFilterLabel(*distMask, *MaskFiltersAnode[p], anode_particle_labels[p], false, -1, CellMode::FULL, Electrode::ANODE);

    //     std::ostringstream name;

    //     name << "MaskFilter_anode_label_"
    //         << std::abs(anode_particle_labels[p])
    //         << ".gf";
    // }

    // // -------------------------------------------------
    // // Cathode particle masks
    // // -------------------------------------------------

    // MaskFiltersCathode.clear();
    // MaskFiltersCathode.resize(cathode_particle_labels.size());

    // for (int p = 0;
    //     p < static_cast<int>(
    //         cathode_particle_labels.size());
    //     ++p)
    // {
    //     MaskFiltersCathode[p] = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    //     ComputePDEFilterLabel(*distMask, *MaskFiltersCathode[p], cathode_particle_labels[p], false, -1, CellMode::FULL, Electrode::CATHODE);

    //     std::ostringstream name;
    //     name << "MaskFilter_cathode_label_"
    //         << cathode_particle_labels[p]
    //         << ".gf";

    //     // MaskFiltersCathode[p]->SaveAsOne(name.str().c_str());
    // }

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

    distMask = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    distMaskSigned = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilterAnode = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilterCathode = std::make_unique<mfem::ParGridFunction>(parfespace.get());
    MaskFilterPse = std::make_unique<mfem::ParGridFunction>(parfespace.get());

    *distMask = 0.0;
    *distMaskSigned = 0.0;
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
    ComputePDEFilter(*distMask, *MaskFilterAnode, 0, CellMode::FULL, Electrode::ANODE);
    ComputePDEFilter(*distMask, *MaskFilterCathode, 0, CellMode::FULL, Electrode::CATHODE);
    ComputePDEFilter(*distMask, *MaskFilterPse, 1, CellMode::FULL, Electrode::ANODE);

    for (std::size_t k = 0; k < anode_particle_labels.size(); ++k)
    {
        ComputePDEFilterLabel(*distMask, *MaskFiltersAnode[k], anode_particle_labels[k], false, -1, CellMode::FULL, Electrode::ANODE);
    }

    for (std::size_t k = 0; k < cathode_particle_labels.size(); ++k)
    {
        ComputePDEFilterLabel(*distMask, *MaskFiltersCathode[k], cathode_particle_labels[k], false, -1, CellMode::FULL, Electrode::CATHODE);
    }
}

void Initialize_Geometry::FullCellAMR()
{
    if (cfg.amr_levels <= 0)
    {
        return;
    }

    mfem::ParGridFunction temporary_dist(parfespace.get());
    mfem::ParGridFunction temporary_anode(parfespace.get());
    mfem::ParGridFunction temporary_cathode(parfespace.get());
    mfem::ParGridFunction temporary_total(parfespace.get());

    temporary_dist = 0.0;
    temporary_anode = 0.0;
    temporary_cathode = 0.0;
    temporary_total = 0.0;

    ComputePDEFilter(temporary_dist, temporary_anode, 0, CellMode::FULL, Electrode::ANODE);
    ComputePDEFilter(temporary_dist, temporary_cathode, 0, CellMode::FULL, Electrode::CATHODE);

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

        temporary_dist.Update();
        temporary_anode.Update();
        temporary_cathode.Update();
        temporary_total.Update();

        ComputePDEFilter(temporary_dist, temporary_anode, 0, CellMode::FULL, Electrode::ANODE);
        ComputePDEFilter(temporary_dist, temporary_cathode, 0, CellMode::FULL, Electrode::CATHODE);

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

    /*
     * A regular integer sort gives:
     *
     * -3, -2, -1
     *
     * Sort by absolute value to get:
     *
     * -1, -2, -3
     */
    std::sort(
        anode_particle_labels.begin(),
        anode_particle_labels.end(),
        [](const int lhs, const int rhs)
        {
            return std::abs(lhs) <
                   std::abs(rhs);
        });

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
            std::cout
                << "[Initialize_Geometry] "
                << "Combining particles separately within "
                << "each full-cell electrode.\n";
        }
    }

    if (myid == 0)
    {
        std::cout
            << "[Initialize_Geometry] "
            << "Anode particle labels: ";

        for (const int label :
             anode_particle_labels)
        {
            std::cout << label << " ";
        }

        std::cout
            << "\n[Initialize_Geometry] "
            << "Cathode particle labels: ";

        for (const int label :
             cathode_particle_labels)
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
    if (voxelData.empty() ||
        voxelData[0].empty() ||
        voxelData[0][0].empty())
    {
        throw std::invalid_argument(
            "InitializeGlobalMesh: voxel data is empty.");
    }

    if (myid == 0)
    {
        std::cout
            << "Creating global mesh from voxel data "
            << "already stored in memory.\n";
    }

    // Copy merged data into the class member.
    tiffData = voxelData;

    globalMesh = CreateGlobalMeshFromTiffData(tiffData);

    globalMesh->EnsureNCMesh(true);

    if (globalMesh->GetNE() <= 0)
    {
        throw std::runtime_error(
            "InitializeGlobalMesh: generated mesh has no elements.");
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
    // const std::string meshFileStr(meshFile);

    // if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) != "tif")
    // {
    //     mfem::mfem_error(
    //         "AssignGlobalValues only supports TIFF files.");
    // }

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << "Reading TIFF file for voxel data\n";
    }

    if (!globalfespace)
    {
        throw std::runtime_error(
            "Global finite element space must be initialized "
            "before assigning global values.");
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

    // Map local to global element indices
    parallelMesh->GetGlobalElementIndices(E_L2G);

    // SetupPinnedDOF(*parfespace);

    gVTX.SetSize(nC);
    VTX.SetSize(nC);

    // Determine file type based on extension
    // std::string meshFileStr(meshFile);  // Convert to std::string
    // if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "tif") {
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

    // SaveTiffDataToPGM(tiffData, "tiff_debug.pgm");

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

        std::cout << "  element size after coarsening = "
                  << cfg.dh * coarsen << "\n";
    }

    return mesh;
}

void Initialize_Geometry::PrintMeshInfo() {
    
    if (!parallelMesh) {
        std::cout << "Parallel mesh not initialized.\n";
        return;
    }

}

void Initialize_Geometry::SaveTiffDataToPGM(
    const std::vector<std::vector<std::vector<int>>> &data,
    const std::string &filename)
{
    if (data.empty() || data[0].empty() || data[0][0].empty()) {
        std::cerr << "SaveTiffDataToPGM: empty data\n";
        return;
    }

    const auto &img = data[0];
    const int height = static_cast<int>(img.size());
    const int width  = static_cast<int>(img[0].size());

    // int max_label = 0;
    // for (const auto &row : img) {
    //     for (int label : row) {
    //         max_label = std::max(max_label, label);
    //     }
    // }

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
                val = static_cast<unsigned char>(
                    std::round(255.0 * (label - minLabel) / (maxLabel - minLabel))
                );
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


void Initialize_Geometry::ComputePDEFilter(mfem::ParGridFunction &dist, mfem::ParGridFunction &filt_gf, int mode, sim::CellMode cell_mode, sim::Electrode electrode)

{
    MFEM_VERIFY(parallelMesh, "parallelMesh is not initialized.");
    MFEM_VERIFY(parfespace, "parfespace is not initialized.");
    MFEM_VERIFY(dist.ParFESpace() == parfespace.get(), "dist must be on parfespace.");
    MFEM_VERIFY(filt_gf.ParFESpace() == parfespace.get(), "filt_gf must be on parfespace.");
    MFEM_VERIFY(parfespace_dg, "parfespace_dg is not initialized.");

    // double dx;
    // dx = parallelMesh->GetElementSize(0); // assuming uniform mesh

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

                    if (mode == 0)
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
                    else if (mode == 1)
                    {
                        fg[idx] = (label == 0) ? 1 : 0;
                    }
                    else
                    {
                        MFEM_ABORT("ComputePDEFilter: mode must be 0 (electrode) or 1 (electrolyte).");
                    }
                }
            }
        }

        if (nz == 1)
        {
            if (mode == 0)
            {
                const int collectorSide = (electrode == Electrode::ANODE) ? 0 : 1;
                KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, collectorSide);
            }
            else
            {
                if (cell_mode == CellMode::FULL)
                {
                    // KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, true, -1);
                    KeepOnlyElectrolyteTouchingBothElectrodes_2D(fg, tiffData, nx, ny, eight_conn);
                }
                else
                {
                    const int electrolyteSide = (electrode == Electrode::ANODE) ? 1 : 0;
                    KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, electrolyteSide);
                    // KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, true, -1);

                }
            }
        }
        else
        {
            if (mode == 0)
            {
                const int collectorFace = (electrode == Electrode::ANODE) ? 0 : 1;
                KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, false, collectorFace);
            }
            else
            {
                if (cell_mode == CellMode::FULL)
                {
                    // KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, true, -1);
                    KeepOnlyElectrolyteTouchingBothElectrodes_3D(fg, tiffData, nx, ny, nz, twenty_six);
                }
                else
                {
                    const int electrolyteFace = (electrode == Electrode::ANODE) ? 1 : 0;
                    KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, false, electrolyteFace);
                }
            }
        }
    }

    // Broadcast full mask to all ranks
    MPI_Bcast(fg.data(), (int)fg.size(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    mfem::ParGridFunction ls_coeff_dg(parfespace_dg.get());
    mfem::ParGridFunction filt_dg(parfespace_dg.get());

    ls_coeff_dg = 0.0;
    filt_dg     = 0.0;

    struct FGCoeffND : public mfem::Coefficient
    {
        int nx, ny, nz;
        int dim;
        double x0,y0,z0, dxp,dyp,dzp;
        const std::vector<uint8_t> *fg;

        FGCoeffND(int nx_, int ny_, int nz_, mfem::ParMesh &pmesh, const std::vector<uint8_t> &fg_)
            : nx(nx_), ny(ny_), nz(nz_), fg(&fg_)
        {
            dim = pmesh.Dimension();
            mfem::Vector bbmin, bbmax;
            pmesh.GetBoundingBox(bbmin, bbmax);

            x0 = bbmin(0); y0 = bbmin(1);
            dxp = (bbmax(0) - bbmin(0)) / (nx - 1);
            dyp = (bbmax(1) - bbmin(1)) / (ny - 1);

            if (dim == 3)
            {
                z0  = bbmin(2);
                dzp = (bbmax(2) - bbmin(2)) / (nz - 1);
            }
            else
            {
                z0 = 0.0; dzp = 1.0;
            }
        }

        double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
        {
            mfem::Vector X;
            T.Transform(ip, X);

            int i = (int)std::floor((X(0) - x0)/dxp + 0.5);
            int j = (int)std::floor((X(1) - y0)/dyp + 0.5);
            i = std::max(0, std::min(nx-1, i));
            j = std::max(0, std::min(ny-1, j));

            int k = 0;
            if (dim == 3)
            {
                k = (int)std::floor((X(2) - z0)/dzp + 0.5);
                k = std::max(0, std::min(nz-1, k));
            }

            const int idx = i + nx*j + nx*ny*k;
            return (*fg)[idx] ? +1.0 : -1.0;
        }
    };

    FGCoeffND fgcoef(nx, ny, nz, *parallelMesh, fg);
    ls_coeff_dg.ProjectCoefficient(fgcoef); 

    // ------------------ PDEFilter ------------------
    const double filter_weight = 3 * cfg.dh;
    mfem::common::PDEFilter filter(*parallelMesh, filter_weight);
    filter.Filter(ls_coeff_dg, filt_dg);

    for (int i = 0; i < filt_dg.Size(); i++)
    {
        filt_dg(i) = 0.5*(filt_dg(i) + 1.0);
    }

    // filt_gf.ProjectGridFunction(filt_dg);

    mfem::GridFunctionCoefficient ls_filt_coeff(&filt_dg);

    filt_gf.ProjectCoefficient(ls_filt_coeff);

    // Explicitly enforce nonconforming H1 constraints.
    mfem::Vector true_values;
    filt_gf.GetTrueDofs(true_values);
    filt_gf.SetFromTrueDofs(true_values);
}

void Initialize_Geometry::ComputePDEFilterLabel(mfem::ParGridFunction &dist, mfem::ParGridFunction &filt_gf, int target_label,
                                                bool keep_boundary_connected, int seed_side_or_face, sim::CellMode cell_mode, sim::Electrode electrode)
{
    MFEM_VERIFY(parallelMesh, "parallelMesh is not initialized.");
    MFEM_VERIFY(parfespace, "parfespace is not initialized.");
    // MFEM_VERIFY(Vox, "Vox is not initialized (need .tif path + MapGlobalToLocal).");
    MFEM_VERIFY(dist.ParFESpace() == parfespace.get(), "dist must be on parfespace.");
    MFEM_VERIFY(filt_gf.ParFESpace() == parfespace.get(), "filt_gf must be on parfespace.");
    MFEM_VERIFY(parfespace_dg, "parfespace_dg is not initialized.");

    // const double dx = parallelMesh->GetElementSize(0);

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
                KeepOnlyConnectedToBoundary_2D(
                    connected_solid,
                    nx,
                    ny,
                    eight_conn,
                    false,
                    seed_side_or_face
                );
            }
            else
            {
                KeepOnlyConnectedToBoundary_3D(
                    connected_solid,
                    nx,
                    ny,
                    nz,
                    twenty_six,
                    false,
                    seed_side_or_face
                );
            }
        }

        // Build this particle-group mask, but only inside the
        // collector-connected solid network.
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

    // if (rank == 0)
    // {
    //     for (int k = 0; k < nz; ++k)
    //     for (int j = 0; j < ny; ++j)
    //     for (int i = 0; i < nx; ++i)
    //     {
    //         const int idx = i + nx*j + nx*ny*k;
    //         // fg[idx] = (tiffData[k][j][i] == target_label) ? 1 : 0;
    //         const int label = tiffData[k][j][i];

    //         if (!cfg.combine_particle_groups)
    //         {
    //             fg[idx] = (label == target_label) ? 1 : 0;
    //         }
    //         else if (cell_mode == CellMode::HALF)
    //         {
    //             fg[idx] = (label > 0) ? 1 : 0;
    //         }
    //         else if (electrode == Electrode::ANODE)
    //         {
    //             fg[idx] = (label < 0) ? 1 : 0;
    //         }
    //         else if (electrode == Electrode::CATHODE)
    //         {
    //             fg[idx] = (label > 0) ? 1 : 0;
    //         }
    //         else
    //         {
    //             MFEM_ABORT("ComputePDEFilterLabel: invalid electrode.");
    //         }
    //     }

    //     if (keep_boundary_connected)
    //     {
    //         if (nz == 1)
    //         {
    //             if (seed_side_or_face < 0)
    //                 KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, true, -1);
    //             else
    //                 KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, seed_side_or_face);
    //         }
    //         else
    //         {
    //             if (seed_side_or_face < 0)
    //                 KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, true, -1);
    //             else
    //                 KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, false, seed_side_or_face);
    //         }
    //     }
    // }

    MPI_Bcast(fg.data(), (int)fg.size(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    mfem::ParGridFunction ls_coeff_dg(parfespace_dg.get());
    mfem::ParGridFunction filt_dg(parfespace_dg.get());

    ls_coeff_dg = 0.0;
    filt_dg = 0.0;

    struct FGCoeffND : public mfem::Coefficient
    {
        int nx, ny, nz;
        int dim;
        double x0, y0, z0, dxp, dyp, dzp;
        const std::vector<uint8_t> *fg;

        FGCoeffND(int nx_, int ny_, int nz_, int dim_,
                  double x0_, double y0_, double z0_,
                  double dxp_, double dyp_, double dzp_,
                  const std::vector<uint8_t> *fg_)
            : nx(nx_), ny(ny_), nz(nz_), dim(dim_),
              x0(x0_), y0(y0_), z0(z0_),
              dxp(dxp_), dyp(dyp_), dzp(dzp_), fg(fg_) {}

        double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
        {
            mfem::Vector x;
            T.Transform(ip, x);

            int i = (int)std::floor((x(0) - x0) / dxp + 0.5);
            int j = (int)std::floor((x(1) - y0) / dyp + 0.5);
            int k = 0;
            if (dim == 3) { k = (int)std::floor((x(2) - z0) / dzp + 0.5); }

            i = std::max(0, std::min(nx-1, i));
            j = std::max(0, std::min(ny-1, j));
            k = std::max(0, std::min(nz-1, k));

            const int idx = i + nx*j + nx*ny*k;
            return ((*fg)[idx] > 0) ? 1.0 : -1.0;
        }
    };

    const int dim = parallelMesh->Dimension();
    mfem::Vector bb_min, bb_max;
    parallelMesh->GetBoundingBox(bb_min, bb_max);

    const double sx = bb_max(0) - bb_min(0);
    const double sy = bb_max(1) - bb_min(1);
    const double sz = (dim == 3) ? (bb_max(2) - bb_min(2)) : cfg.dh;

    const double x0 = bb_min(0);
    const double y0 = bb_min(1);
    const double z0 = (dim == 3) ? bb_min(2) : 0.0;

    const double dxp = (nx > 1) ? sx / (nx - 1) : cfg.dh;
    const double dyp = (ny > 1) ? sy / (ny - 1) : cfg.dh;
    const double dzp = (dim == 3 && nz > 1) ? sz / (nz - 1) : cfg.dh;

    FGCoeffND fg_coeff(nx, ny, nz, dim, x0, y0, z0, dxp, dyp, dzp, &fg);
    ls_coeff_dg.ProjectCoefficient(fg_coeff);

    // ------------------ PDEFilter ------------------
    const double filter_weight = 3 * cfg.dh;
    mfem::common::PDEFilter filter(*parallelMesh, filter_weight);
    filter.Filter(ls_coeff_dg, filt_dg);

    for (int i = 0; i < filt_dg.Size(); i++)
    {
        filt_dg(i) = 0.5*(filt_dg(i) + 1.0);
    }

    // filt_gf.ProjectGridFunction(filt_dg);

    mfem::GridFunctionCoefficient ls_filt_coeff(&filt_dg);

    filt_gf.ProjectCoefficient(ls_filt_coeff);

    // Explicitly enforce nonconforming H1 constraints.
    mfem::Vector true_values;
    filt_gf.GetTrueDofs(true_values);
    filt_gf.SetFromTrueDofs(true_values);
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
