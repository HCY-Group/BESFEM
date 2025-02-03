#include "Initialize_Geometry.hpp"
#include "Constants.hpp"
#include "readtiff.h"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

// Constructor
Initialize_Geometry::Initialize_Geometry() 
{}

// Destructor
Initialize_Geometry::~Initialize_Geometry() {}

void Initialize_Geometry::InitializeMesh(const char* meshFile, MPI_Comm comm, int order) {

    // Initialize the global mesh
    InitializeGlobalMesh(Constants::mesh_file);

    // Initialize the parallel mesh
    InitializeParallelMesh(MPI_COMM_WORLD);

    // Set up the finite element space
    SetupFiniteElementSpace(Constants::order);

    // Set up the parallel finite element space
    SetupParFiniteElementSpace(Constants::order);

    // Assign the global values
    AssignGlobalValues(Constants::mesh_file);

    // Map the global values to the local
    MapGlobalToLocal(Constants::mesh_file);

    // Print out information relative to the mesh
    PrintMeshInfo();

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
        std::cout << "Creating global mesh using .mesh file" << std::endl;
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
    auto fec = new mfem::H1_FECollection(order, globalMesh->Dimension());
    globalfespace = std::make_shared<mfem::FiniteElementSpace>(globalMesh.get(), fec);
}

// Function to set up finite element space on parallel mesh
void Initialize_Geometry::SetupParFiniteElementSpace(int order) {
    if (!parallelMesh) {
        throw std::runtime_error("Parallel mesh must be initialized before setting up FE space.");
    }
    auto fec = new mfem::H1_FECollection(order, parallelMesh->Dimension());
    parfespace = std::make_shared<mfem::ParFiniteElementSpace>(parallelMesh.get(), fec);
}


void Initialize_Geometry::AssignGlobalValues(const char* meshFile) {
    std::string meshFileStr(meshFile);  // Convert to std::string
    if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "tif") {
        // Process the .tif file
        cout << "Reading .tif file for voxel data" << endl;

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

    } else if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "mesh") {
        // Process the .mesh file
        cout << "Reading .mesh file for global distance function" << endl;
        ifstream myfile(meshFile);
        if (myfile.is_open()) {
            // Assuming that the mesh data has a similar structure
            gDsF = make_unique<mfem::GridFunction>(parfespace.get());
            Onm = gDsF->Size();
            for (int gi = 0; gi < Onm; gi++) {
                myfile >> (*gDsF)(gi);
            }
            myfile.close();
        } else {
            cerr << "Failed to open .mesh file" << endl;
        }
    } else {
        cerr << "Unsupported file type" << endl;
    }
}

void Initialize_Geometry::MapGlobalToLocal(const char* meshFile) {
    
    if (!parallelMesh) {
        throw std::runtime_error("Parallel mesh must be initialized before calculating element volumes.");
    }

    if (!globalMesh) {
        throw std::runtime_error("Global mesh must be initialized before setting up FE space.");
    }
        
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
        cout << "Reading .tif file for mapping global to local grid function" << endl;

        Vox = std::make_unique<mfem::ParGridFunction>(parfespace.get());

        // Iterate over elements and map global to local
        for (ei = 0; ei < nE; ei++) {
            gei = E_L2G[ei];

            globalMesh->GetElementVertices(gei, gVTX);
            parallelMesh->GetElementVertices(ei, VTX);

            for (int vi = 0; vi < nC; vi++) {
                (*this->Vox)(VTX[vi]) = (*this->gVox)(gVTX[vi]);
            }
        }

    } else if (meshFileStr.substr(meshFileStr.find_last_of(".") + 1) == "mesh") {
        // Handle .mesh file
        cout << "Reading .mesh file for mapping global to local grid function" << endl;

        // Read global distance function
        gDsF = make_unique<mfem::GridFunction>(parfespace.get());
        Onm = gDsF->Size();
        ifstream myfile(Constants::dsF_file);
        for (int gi = 0; gi < Onm; gi++) {
            myfile >> (*gDsF)(gi);
        }
        myfile.close();

        // Assuming dsF is a ParGridFunction for the mesh file
        dsF = make_unique<mfem::ParGridFunction>(parfespace.get());

        // Map local distance function from global one
        for (ei = 0; ei < nE; ei++) {
            gei = E_L2G[ei];

            globalMesh->GetElementVertices(gei, gVTX);
            parallelMesh->GetElementVertices(ei, VTX);

            for (int vi = 0; vi < nC; vi++) {
                (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]);
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
	
	cout << "tiff size: "            << tiffData.size()         << endl;	
	cout << "tiff[0] size: "         << tiffData[0].size()      << endl;
	cout << "tiff[0][0] size: "      << tiffData[0][0].size()   << endl;
	cout << "tiffdata[0][0][0] = "   << tiffData[0][0][0]       << endl;

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

    std::cout << "Number of vertices: " << nV << "\n";
    std::cout << "Number of elements: " << nE << "\n";
}

// void Mesh::OutputToParaview(const std::string &fileName, const std::string &varName, mfem::GridFunction *gf) {
//     auto pd = new mfem::ParaViewDataCollection(fileName, gf->FESpace()->GetMesh());
//     pd->RegisterField(varName, gf);
//     pd->SetLevelsOfDetail(1);
//     pd->SetDataFormat(mfem::VTKFormat::BINARY);
//     pd->SetHighOrderOutput(true);
//     pd->Save();
//     delete pd;
// }