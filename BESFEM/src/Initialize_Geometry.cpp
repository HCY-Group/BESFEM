#include "../include/Initialize_Geometry.hpp"
#include "../inputs/Constants.hpp"
#include "../include/readtiff.h"
#include "../include/SimTypes.hpp"
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

// Half Cell
void Initialize_Geometry::InitializeMesh(const char* meshFile, const char* distanceFile, MPI_Comm comm, int order) {

    // Adjust distance file
    AdjustDistanceFile(distanceFile);
    
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

    // Print out information relative to the mesh
    PrintMeshInfo();

}

// Full Cell
void Initialize_Geometry::InitializeMesh(const char* meshFile, const char* distanceFileA, const char* distanceFileC, MPI_Comm comm, int order) {

    // Adjust distance file
    AdjustDistanceFile(distanceFileA); // for anode
    AdjustDistanceFile(distanceFileC); // for cathode

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

void Initialize_Geometry::AdjustDistanceFile(const char* distanceFile)
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

            std::cout << "[AdjustDistanceFile] Before: min=" << *min_it
                      << " max=" << *max_it << " max|v|=" << max_abs << "\n";

            if (max_abs > 1.0) {
                // 1) write backup with original data
                std::string backup = std::string(distanceFile) + ".orig";
                {
                    std::ofstream bout(backup, std::ios::trunc);
                    if (!bout) {
                        throw std::runtime_error("Failed to create backup file: " + backup);
                    }
                    bout << std::setprecision(10);
                    for (double x : values) bout << x << '\n';
                }
                std::cout << "[AdjustDistanceFile] Wrote original data to backup: " << backup << "\n";

                // 2) scale and overwrite original file
                std::cout << "[AdjustDistanceFile] Scaling values by dh=" << std::setprecision(10) << Constants::dh << " and overwriting "
                          << distanceFile << " ...\n";
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
                std::cout << "[AdjustDistanceFile] After: min=" << *min2
                          << " max=" << *max2 << " max|v|=" << max_abs2 << "\n";
            } else {
                std::cout << "[AdjustDistanceFile] No scaling needed (all |v| <= 1). File unchanged.\n";
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
        std::cout << "Creating global mesh using .tif file" << std::endl;
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
}

// Function to initialize the parallel mesh
void Initialize_Geometry::InitializeParallelMesh(MPI_Comm comm) {
    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before creating a parallel mesh.");
    }
    parallelMesh = std::make_unique<mfem::ParMesh>(comm, *globalMesh);
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
                    cerr << "Error: Distance file has fewer than four header lines" << endl;
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
    mfem::Array<HYPRE_BigInt> E_L2G;
    parallelMesh->GetGlobalElementIndices(E_L2G);

    mfem::Array<int> gVTX(nC);    // global indices of corner vertices
    mfem::Array<int> VTX(nC);     // local indices of corner vertices

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

	cout << "reading tiff file" << endl;
	Constraints args;
	//TODO: The code works with serial 2d, parallel 2d, and serial 3d, but not parallel 3d
	args.Depth_begin = 0;	//only read in one slice for 2D data
	args.Depth_end = 1;	//only read in one slice for 2D data
	// get a smaller subset so it runs faster
	args.Row_begin    = 0;
	args.Row_end      = 80;
	args.Column_begin = 0;
	args.Column_end   = 120;
	TIFFReader reader(meshFile,args);
	reader.readinfo();
	std::vector<std::vector<std::vector<int>>> tiffData;
	tiffData = reader.getImageData();

    return tiffData;
}

// Create a global MFEM mesh from voxel data extracted from .tif file
std::unique_ptr<mfem::Mesh> Initialize_Geometry::CreateGlobalMeshFromTiffData(const std::vector<std::vector<std::vector<int>>>& tiffData) {
    int nz = tiffData.size(); // depth dimension
    int ny = tiffData[0].size(); // row dimension
    int nx = tiffData[0][0].size(); // column dimension
    double sx = nx;  // make dx = 1 // size in x direction
    double sy = ny;  // make dy = 1 // size in y direction
    double sz = nz;  // make dz = 1 // size in z direction
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
}

void Initialize_Geometry::PrintMeshInfo() {
    
    if (!parallelMesh) {
        std::cout << "Parallel mesh not initialized.\n";
        return;
    }

}

void Initialize_Geometry::SetupBoundaryConditions(CellMode mode, Electrode electrode) {

    int dim = parallelMesh->Dimension();

    if (dim == 3) {

    std::cout << "Setting up boundary conditions for 3D mesh" << std::endl;

    if (mode == CellMode::HALF && electrode == Electrode::ANODE) {
        std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;

        // East Neumann Boundary Condition
        nbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_e_bdr = 0;
        nbc_e_bdr[2] = 1; 

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());

        nbc_bdr = nbc_e_bdr;
        dbc_bdr = dbc_e_bdr;


    } else if (mode == CellMode::HALF && electrode == Electrode::CATHODE) {
        std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;

        // West Neumann Boundary Condition
        nbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_w_bdr = 0;
        nbc_w_bdr[0] = 1; 

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        // East Dirichlet Boundary Condition 
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr = nbc_w_bdr;
        dbc_bdr = dbc_w_bdr;

    } else {
        std::cout << "Setting up boundary conditions for Full Cell" << std::endl;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_bdr = 0;

        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr = 0;

    }
    
    // // DISK
    // // Boundary attributes for Neumann BC on the west boundary for CnE
    // nbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
    // nbc_w_bdr = 0;
    // nbc_w_bdr[2] = 1; // east

    // // Dirichlet BC on the east boundary for CnP & phP
    // dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
    // dbc_e_bdr = 0; 
    // dbc_e_bdr[0] = 1;  // Applying Dirichlet BC to the west boundary

    // // Extract essential true DOFs (Dirichlet BCs) on the east boundary
    // ess_tdof_list_e.SetSize(0);
    // parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

    // // Dirichlet BC on the west boundary for CnE & phE
    // dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
    // dbc_w_bdr = 0;
    // dbc_w_bdr[2] = 1;   // east

    // // Extract essential true DOFs (Dirichlet BCs) on the west boundary
    // ess_tdof_list_w.SetSize(0); 
    // parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

    } else if (dim == 2) {

    std::cout << "Setting up boundary conditions for 2D mesh" << std::endl;

    if (mode == CellMode::HALF && electrode == Electrode::ANODE) {
        std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;

        // East Neumann Boundary Condition
        nbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_e_bdr = 0;
        nbc_e_bdr[2] = 1; 

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());

        nbc_bdr = nbc_e_bdr;
        dbc_bdr = dbc_e_bdr;


    } else if (mode == CellMode::HALF && electrode == Electrode::CATHODE) {
        std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;

        // West Neumann Boundary Condition
        nbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_w_bdr = 0;
        nbc_w_bdr[0] = 1; 

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        // East Dirichlet Boundary Condition 
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr = nbc_w_bdr;
        dbc_bdr = dbc_w_bdr;

    } else {
        std::cout << "Setting up boundary conditions for Full Cell" << std::endl;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_bdr = 0;

        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr = 0;

    }


}

}

