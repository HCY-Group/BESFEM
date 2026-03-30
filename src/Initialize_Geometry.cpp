#include "../include/Initialize_Geometry.hpp"
#include "../include/Constants.hpp"
#include "../include/readtiff.h"
#include "../include/SimTypes.hpp"
#include "../include/dist_solver.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;
using sim::CellMode;
using sim::Electrode;


// Constructor
Initialize_Geometry::Initialize_Geometry() 
{}

// Destructor
Initialize_Geometry::~Initialize_Geometry() {}


#include <queue>
#include <cstdint>

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



// Half Cell
void Initialize_Geometry::InitializeMesh(const char* meshFile, const char* distanceFile, const char* mesh_type, MPI_Comm comm, int order) {

    myid = mfem::Mpi::WorldRank();

    // Adjust distance file
    AdjustDistanceFile(distanceFile, mesh_type);
    
    // Initialize the global mesh
    InitializeGlobalMesh(meshFile);

    // Initialize the parallel mesh
    InitializeParallelMesh(MPI_COMM_WORLD);

    // Set up the finite element space
    SetupFiniteElementSpace(order);

    // Set up the parallel finite element space
    SetupParFiniteElementSpace(order);

    // Assign the global values
    AssignGlobalValues(meshFile, distanceFile, gDsF);

    // Map the global values to the local
    MapGlobalToLocal(meshFile);

    // distance for tiff files
    std::string meshFileStr(meshFile);
    if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "tif")
    {
        distMask      = std::make_unique<mfem::ParGridFunction>(parfespace.get());
        distMaskSigned= std::make_unique<mfem::ParGridFunction>(parfespace.get());
        MaskFilter   = std::make_unique<mfem::ParGridFunction>(parfespace.get());
        MaskFilterPse = std::make_unique<mfem::ParGridFunction>(parfespace.get());

        const int solver_type = 0;
        const double t_param = 1.0;

        ComputePDEFilter(*distMask, *MaskFilter, /*mode=*/0); 
        ComputePDEFilter(*distMask, *MaskFilterPse, /*mode=*/1);

        if (mfem::Mpi::WorldRank() == 0) { std::cout << "ComputePDEFilter done.\n"; }

        MaskFilter->SaveAsOne("MaskFilter.gf");
        MaskFilterPse->SaveAsOne("MaskFilter_pse.gf");

    }

    // Print out information relative to the mesh
    PrintMeshInfo();

    globalMesh->Save("gmesh");


}

// Full Cell
void Initialize_Geometry::InitializeMesh(const char* meshFile, const char* distanceFileA, const char* distanceFileC, const char* mesh_type, MPI_Comm comm, int order) {

    myid = mfem::Mpi::WorldRank();

    // Adjust distance file
    AdjustDistanceFile(distanceFileA, mesh_type); // for anode
    AdjustDistanceFile(distanceFileC, mesh_type); // for cathode

    // Initialize the global mesh
    InitializeGlobalMesh(meshFile);

    // Initialize the parallel mesh
    InitializeParallelMesh(MPI_COMM_WORLD);

    // Set up the finite element space
    SetupFiniteElementSpace(order);

    // Set up the parallel finite element space
    SetupParFiniteElementSpace(order);

    // Assign the global values
    AssignGlobalValues(meshFile, distanceFileA, gDsF_A); // for anode
    AssignGlobalValues(meshFile, distanceFileC, gDsF_C); // for cathode

    // Map the global values to the local
    MapGlobalToLocal(meshFile);

    // Print out information relative to the mesh
    PrintMeshInfo();

}

void Initialize_Geometry::AdjustDistanceFile(const char* distanceFile, const char* mesh_type)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        std::ifstream in(distanceFile);
        if (!in) {
            throw std::runtime_error("Could not open distance file: " + std::string(distanceFile));
        }

        std::vector<double> values;
        values.reserve(1<<20); // pre-alloc for big files (optional)

        double v;
        while (in >> v) values.push_back(v);
        in.close();

        if (values.empty()) {
            std::cerr << "[AdjustDistanceFile] No values read from " << distanceFile << " — leaving unchanged.\n";
        } else {
            auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
            double max_abs = std::max(std::abs(*min_it), std::abs(*max_it));

            // std::cout << "[AdjustDistanceFile] Before: min=" << *min_it
            //           << " max=" << *max_it << " max|v|=" << max_abs << "\n";

            if (strcmp(mesh_type, "ml") == 0 && max_abs > 1.0) {
                // write backup with original data
                std::string backup = std::string(distanceFile) + ".orig";
                {
                    std::ofstream bout(backup, std::ios::trunc);
                    if (!bout) {
                        throw std::runtime_error("Failed to create backup file: " + backup);
                    }
                    bout << std::setprecision(10);
                    for (double x : values) bout << x << '\n';
                }

                if (mfem::Mpi::WorldRank() == 0) { std::cout << "[AdjustDistanceFile] Wrote original data to backup: " << backup << "\n"; }

                // scale and overwrite original file
                if (mfem::Mpi::WorldRank() == 0) { std::cout << "[AdjustDistanceFile] Scaling values by dh=" << std::setprecision(10) << Constants::dh << " and overwriting "
                          << distanceFile << " ...\n"; }
                for (double &x : values) x *= Constants::dh;

                std::ofstream out(distanceFile, std::ios::trunc);
                if (!out) {
                    throw std::runtime_error("Failed to open distance file for overwrite: " + std::string(distanceFile));
                }
                out << std::setprecision(10);
                for (double x : values) out << x << '\n';
                out.close();

                // quick preview after
                auto [min2, max2] = std::minmax_element(values.begin(), values.end());
                double max_abs2 = std::max(std::abs(*min2), std::abs(*max2));
                // std::cout << "[AdjustDistanceFile] After: min=" << *min2
                //           << " max=" << *max2 << " max|v|=" << max_abs2 << "\n";
            } else {
                if (mfem::Mpi::WorldRank() == 0) { std::cout << "[AdjustDistanceFile] No scaling needed (all |v| <= 1 or voxel). File unchanged.\n"; }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
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
    else if (fileExtension == "mesh") {
        if (mfem::Mpi::WorldRank() == 0) // only print on rank 0
        {std::cout << "Creating global mesh using .mesh file" << std::endl;}
        
        globalMesh = std::make_unique<mfem::Mesh>(meshFile);
    } 
    else {
        throw std::invalid_argument("Unsupported file format. Only .tif and .mesh are allowed.");
    }

    // ensure mesh supports non-conforming elements for adaptive refinement
    globalMesh->EnsureNCMesh(true);


    int e = 0;
    mfem::Array<int> vert_ids;
    globalMesh->GetElementVertices(e, vert_ids);

    mfem::Vector v0(globalMesh->GetVertex(vert_ids[0]), globalMesh->SpaceDimension());
    mfem::Vector v1(globalMesh->GetVertex(vert_ids[1]), globalMesh->SpaceDimension());


    double dh1 = v0.DistanceTo(v1);
    if (mfem::Mpi::WorldRank() == 0) { std::cout << "Element size dh = " << dh1 << std::endl;}

    MPI_Barrier(MPI_COMM_WORLD);

}

// Function to initialize the parallel mesh
void Initialize_Geometry::InitializeParallelMesh(MPI_Comm comm) {
    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before creating a parallel mesh.");
    }
    parallelMesh = std::make_shared<mfem::ParMesh>(comm, *globalMesh);
    parallelMesh->SaveAsOne("pmesh");



    // std::cout << "Rank " << myid << " owns "
    //           << parallelMesh->GetNE() << " elements, "
    //           << parallelMesh->GetNV() << " vertices.\n";
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


void Initialize_Geometry::AssignGlobalValues(const char* meshFile, const char* distanceFile, std::unique_ptr<mfem::GridFunction>& gDsF_out) {
    std::string meshFileStr(meshFile);  // Convert to std::string
    
    if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "tif") {
        // Process the .tif file
        if (mfem::Mpi::WorldRank() == 0){ // only print on rank 0 
        cout << "Reading .tif file for voxel data" << endl;
        }

    if (!globalfespace) {
        throw std::runtime_error("Global finite element space (globalfespace) must be initialized before assigning global values.");
    }
        
    gVox = std::make_unique<mfem::GridFunction>(globalfespace.get());

        int nz = tiffData.size();
        int ny = tiffData[0].size();
        int nx = tiffData[0][0].size();
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    int idx = i + nx * j + nx * ny * k;
                    (*this->gVox)[idx] = tiffData[k][j][i];
                }
            }
        }
    if (mfem::Mpi::WorldRank() == 0) { // only print on rank 0
    cout << "Reading .dsF file for global distance function for tif case" << endl;
    }
        gDsF_out = make_unique<mfem::GridFunction>(globalfespace.get());
        std::ifstream myfile(distanceFile);
        if (myfile.is_open()) {
            // Skip the first four lines
            string line;
            for (int i = 0; i < 4; i++) {
                if (!getline(myfile, line)) {
                    if (mfem::Mpi::WorldRank() == 0) {cerr << "Warning: Distance file has fewer than four header lines" << endl;}
                    myfile.close();
                    return;
                }
            }
            gDsF_out->Load(myfile, gDsF_out->Size());
            myfile.close();
        } else {
            cerr << "Failed to open distance file" << endl;
        }
        

    } else if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "mesh") {
    
    if (mfem::Mpi::WorldRank() == 0) // only print on rank 0
    { cout << "Reading .dsF file for global distance function for mesh case" << endl; }

        gDsF_out = make_unique<mfem::GridFunction>(globalfespace.get());
        ifstream myfile(distanceFile);
        if (myfile.is_open()) {
            gDsF_out->Load(myfile, gDsF_out->Size());
            myfile.close();
        } else {
            cerr << "Failed to open distance file" << endl;
        }
    }
}   

void Initialize_Geometry::MapGlobalToLocal(const char* meshFile) {
    
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
    std::string meshFileStr(meshFile);  // Convert to std::string
    if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "tif") {
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
            
            if (gDsF) {
                if (!dsF) dsF = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi = 0; vi < nC; ++vi) { (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]); }
            }

            if (gDsF_A) {
                if (!dsF_A) dsF_A = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi=0; vi<nC; ++vi) (*dsF_A)(VTX[vi]) = (*gDsF_A)(gVTX[vi]);
            }
            if (gDsF_C) {
                if (!dsF_C) dsF_C = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi=0; vi<nC; ++vi) (*dsF_C)(VTX[vi]) = (*gDsF_C)(gVTX[vi]);
            }

            if (!gDsF_C && gDsF_A) {
                if (!dsF) dsF = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi=0; vi<nC; ++vi) (*dsF)(VTX[vi]) = (*gDsF_A)(gVTX[vi]);
            }
        }

    } else if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "mesh") {
        // Handle .mesh file
        if (mfem::Mpi::WorldRank() == 0) // only print on rank 0
        {cout << "Reading .mesh file for mapping global to local grid function" << endl;}

        // Map local distance function from global one
        for (ei = 0; ei < nE; ei++) {
            gei = E_L2G[ei];

            globalMesh->GetElementVertices(gei, gVTX);
            parallelMesh->GetElementVertices(ei, VTX);

            if (gDsF) {
                if (!dsF) dsF = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi = 0; vi < nC; ++vi) { (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]); }
            }

            if (gDsF_A) {
                if (!dsF_A) dsF_A = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi=0; vi<nC; ++vi) (*dsF_A)(VTX[vi]) = (*gDsF_A)(gVTX[vi]);
            }
            if (gDsF_C) {
                if (!dsF_C) dsF_C = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi=0; vi<nC; ++vi) (*dsF_C)(VTX[vi]) = (*gDsF_C)(gVTX[vi]);
            }

            if (!gDsF_C && gDsF_A) {
                if (!dsF) dsF = std::make_unique<mfem::ParGridFunction>(parfespace.get());
                for (int vi=0; vi<nC; ++vi) (*dsF)(VTX[vi]) = (*gDsF_A)(gVTX[vi]);
            }

        }
    } else {
        cerr << "Unsupported file type for MapGlobalToLocal" << endl;
    }

}

// Reading .tif file and returning voxel data
std::vector<std::vector<std::vector<int>>> Initialize_Geometry::ReadTiffFile(const char* meshFile) {

	if (mfem::Mpi::WorldRank() == 0) { std::cout << "reading tiff file" << std::endl; }
	Constraints args;
	//TODO: The code works with serial 2d, parallel 2d, and serial 3d, but not parallel 3d
	args.Depth_begin = 0;	//only read in one slice for 2D data
	args.Depth_end = 1;	//only read in one slice for 2D data
	// get a smaller subset so it runs faster
	args.Row_begin    = 0;
	args.Row_end      = 100;
	args.Column_begin = 0;
	args.Column_end   = 100;
	TIFFReader reader(meshFile,args);
	reader.readinfo();
	std::vector<std::vector<std::vector<int>>> tiffData;
	tiffData = reader.getImageData();

    SaveTiffDataToPGM(tiffData, "tiff_debug.pgm");

    return tiffData;
}

// Create a global MFEM mesh from voxel data extracted from .tif file
std::unique_ptr<mfem::Mesh> Initialize_Geometry::CreateGlobalMeshFromTiffData(const std::vector<std::vector<std::vector<int>>>& tiffData) {
    int nz = tiffData.size(); // depth dimension
    int ny = tiffData[0].size(); // row dimension
    int nx = tiffData[0][0].size(); // column dimension
    
    // double sx = nx;  // make dx = 1 // size in x direction
    // double sy = ny;  // make dy = 1 // size in y direction
    // double sz = nz;  // make dz = 1 // size in z direction

    double scale = 2.0e-5;

    double sx = nx * scale;  // make dx = 1 // size in x direction
    double sy = ny * scale;  // make dy = 1 // size in y direction
    double sz = nz * scale;  // make dz = 1 // size in z direction


    bool generate_edges = false; 
    bool sfc_ordering = false; 

    std::unique_ptr<mfem::Mesh> mesh;

    if (nz == 1) {
        mesh = std::make_unique<mfem::Mesh>(
            mfem::Mesh::MakeCartesian2D(nx - 1, ny - 1, mfem::Element::QUADRILATERAL, generate_edges, sx, sy, sfc_ordering)
        );
    } else {
        mesh = std::make_unique<mfem::Mesh>(
            mfem::Mesh::MakeCartesian3D(nx - 1, ny - 1, nz - 1, mfem::Element::HEXAHEDRON, sx, sy, sz, sfc_ordering)
        );
    }

    return mesh;

    // std::cout << "[CreateGlobalMeshFromTiffData] nx=" << nx
    //       << " ny=" << ny << " nz=" << nz << "\n"
    //       << " sx=" << sx << " sy=" << sy << " sz=" << sz << std::endl;

}

void Initialize_Geometry::PrintMeshInfo() {
    
    if (!parallelMesh) {
        std::cout << "Parallel mesh not initialized.\n";
        return;
    }

}

void Initialize_Geometry::SaveTiffDataToPGM(const std::vector<std::vector<std::vector<int>>> &data,
                              const std::string &filename)
{
    if (data.empty() || data[0].empty() || data[0][0].empty()) {
        std::cerr << "SaveTiffDataToPGM: empty data\n";
        return;
    }

    const auto &img = data[0];              // first slice only
    const int height = (int)img.size();     // rows
    const int width  = (int)img[0].size();  // columns

    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Could not open file for writing: " << filename << "\n";
        return;
    }

    // PGM header
    out << "P5\n" << width << " " << height << "\n255\n";

    // Write binary 0 or 255 only
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            unsigned char val;

            if (img[j][i] <= 0)       val = 0;     // black
            else                      val = 255;   // white

            out.write(reinterpret_cast<char*>(&val), 1);
        }
    }

    out.close();
    if (mfem::Mpi::WorldRank() == 0) {std::cout << "Saved binary PGM (0/255) to " << filename << "\n";}
}


void Initialize_Geometry::ComputePDEFilter(mfem::ParGridFunction &dist, mfem::ParGridFunction &filt_gf, int mode)

{
    MFEM_VERIFY(parallelMesh, "parallelMesh is not initialized.");
    MFEM_VERIFY(parfespace, "parfespace is not initialized.");
    MFEM_VERIFY(Vox, "Vox is not initialized (need .tif path + MapGlobalToLocal).");
    MFEM_VERIFY(dist.ParFESpace() == parfespace.get(), "dist must be on parfespace.");
    MFEM_VERIFY(filt_gf.ParFESpace() == parfespace.get(), "filt_gf must be on parfespace.");
    MFEM_VERIFY(parfespace_dg, "parfespace_dg is not initialized.");
    // MFEM_VERIFY(parallelMesh->Dimension() == 2, "This 2D connectivity helper assumes a 2D TIFF slice.");

    double dx;
    dx = parallelMesh->GetElementSize(0); // assuming uniform mesh

    // new psi method?? PDEFilter - Poisson smoothing field

    // // // like doughnut and cheese (1 inside, -1 outside)
    // mfem::ParGridFunction ls_coeff(parfespace.get());
    // for (int i = 0; i < ls_coeff.Size(); i++)
    // {
    //     const double m = (*Vox)(i);            
    // //     // ls_coeff(i) = (m > 0.5) ? 0.0 : +1.0; // define what is inside vs outside (other tif with black outline)
    //     // ls_coeff(i) = (m > 0.5) ? +1.0 : 0.0; // define what is inside vs outside (microstructure tif) HEAT
    //     ls_coeff(i) = (m > 0.5) ? +1.0 : -1.0; // define what is inside vs outside (microstructure tif) P LAP
        

    // }

    MFEM_VERIFY(parallelMesh->Dimension() == 2 || parallelMesh->Dimension() == 3,
            "ComputePDEFilter: mesh must be 2D or 3D.");

    // TIFF sizes
    const int nz = (int)tiffData.size();
    const int ny = (int)tiffData[0].size();
    const int nx = (int)tiffData[0][0].size();

    const int rank = mfem::Mpi::WorldRank();
    const bool eight_conn = false;        // 2D: 4/8
    const bool twenty_six = false;        // 3D: 6/26  (set true if desired)

    // IMPORTANT: fg must cover the whole volume for 3D
    std::vector<uint8_t> fg(nx*ny*nz, 0);

    if (rank == 0)
    {
        // base solid from TIFF: 1 = white/solid
        std::vector<uint8_t> solid_base(nx*ny*nz, 0);
        for (int k=0; k<nz; ++k)
        for (int j=0; j<ny; ++j)
        for (int i=0; i<nx; ++i)
        {
            const int idx = i + nx*j + nx*ny*k;
            solid_base[idx] = (tiffData[k][j][i] > 0) ? 1 : 0;
        }

        if (mode == 0)
        {
            // PSI: solid phase, keep only stuff connected to a boundary
            fg = solid_base;

            if (nz == 1)
            {
                // right boundary in 2D (seed_side=1)
                KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, 1); // right boundary in 2D (seed_side=1)
            }
            else
            {
                // xmax face in 3D (seed_face=1)
                KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, false, 1);
            }
        }
        else if (mode == 1)
        {
            // PSE: electrolyte = NOT solid
            for (int idx=0; idx<nx*ny*nz; ++idx) fg[idx] = solid_base[idx] ? 0 : 1;

            if (nz == 1)
            {
                KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, 0); // left boundary in 2D (seed_side=0)
                // KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, true, -1); // all boundaries
            }
            else
            {
                KeepOnlyConnectedToBoundary_3D(fg, nx, ny, nz, twenty_six, false, 0);
            }
        }
        else
        {
            MFEM_ABORT("ComputePDEFilter: mode must be 0 (psi) or 1 (pse).");
        }
    }

    // Broadcast full mask to all ranks
    MPI_Bcast(fg.data(), (int)fg.size(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);


    // const int nv_loc = parallelMesh->GetNV();
    // const int ny = (int)tiffData[0].size();
    // const int nx = (int)tiffData[0][0].size();
    // const int rank = mfem::Mpi::WorldRank();
    // const bool eight_conn = false;

    mfem::ParGridFunction ls_coeff_dg(parfespace_dg.get());
    mfem::ParGridFunction filt_dg(parfespace_dg.get());

    ls_coeff_dg = 0.0;
    filt_dg     = 0.0;

    // std::vector<uint8_t> fg(nx*ny, 0);

    // if (rank == 0)
    // {
    //     // base solid from TIFF: 1 = white/solid 
    //     std::vector<uint8_t> solid_base(nx*ny, 0);
    //     for (int j=0; j<ny; ++j)
    //     for (int i=0; i<nx; ++i)
    //     {
    //         const int k = i + nx*j;
    //         solid_base[k] = (tiffData[0][j][i] > 0) ? 1 : 0;
    //     }

    //     if (mode == 0)
    //     {
    //         // --- PSI
    //         fg = solid_base;
    //         KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, false, 1); // right
    //     }
    //     else if (mode == 1)
    //     {
    //         // --- PSE
    //         for (int k=0; k<nx*ny; ++k) fg[k] = solid_base[k] ? 0 : 1; // fg=1 for electrolyte 
    //         KeepOnlyConnectedToBoundary_2D(fg, nx, ny, eight_conn, true, -1); // left
    //     }
    //     else
    //     {
    //         MFEM_ABORT("ComputeDistanceFromTiffMask: mode must be 0 (psi) or 1 (pse).");
    //     }
    // }

    // MPI_Bcast(fg.data(), nx*ny, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

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


    // struct FGCoeff : public mfem::Coefficient
    // {
    //     int nx, ny;
    //     double x0, y0, dxp, dyp;
    //     const std::vector<uint8_t> *fg;

    //     FGCoeff(int nx_, int ny_, mfem::ParMesh &pmesh, const std::vector<uint8_t> &fg_)
    //         : nx(nx_), ny(ny_), fg(&fg_)
    //     {
    //         mfem::Vector bbmin, bbmax;
    //         pmesh.GetBoundingBox(bbmin, bbmax);
    //         x0  = bbmin(0);
    //         y0  = bbmin(1);
    //         dxp = (bbmax(0) - bbmin(0)) / (nx - 1);
    //         dyp = (bbmax(1) - bbmin(1)) / (ny - 1);
    //     }

    //     double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
    //     {
    //         mfem::Vector X;
    //         T.Transform(ip, X);

    //         int i = (int)std::floor((X(0) - x0)/dxp + 0.5);
    //         int j = (int)std::floor((X(1) - y0)/dyp + 0.5);

    //         if (i < 0) i = 0; if (i > nx-1) i = nx-1;
    //         if (j < 0) j = 0; if (j > ny-1) j = ny-1;

    //         return (*fg)[i + nx*j] ? +1.0 : -1.0;
    //     }
    // };


    // FGCoeff fgcoef(nx, ny, *parallelMesh, fg);
    // ls_coeff_dg.ProjectCoefficient(fgcoef);

    FGCoeffND fgcoef(nx, ny, nz, *parallelMesh, fg);
    ls_coeff_dg.ProjectCoefficient(fgcoef); 



    // ------------------ PDEFilter ------------------
    const double filter_weight = 3 * dx;
    mfem::common::PDEFilter filter(*parallelMesh, filter_weight);
    filter.Filter(ls_coeff_dg, filt_dg);

    for (int i = 0; i < filt_dg.Size(); i++)
    {
        filt_dg(i) = 0.5*(filt_dg(i) + 1.0);
    }

    mfem::GridFunctionCoefficient ls_filt_coeff(&filt_dg);

    filt_gf.ProjectGridFunction(filt_dg);
}
