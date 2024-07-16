// #include "mfem.hpp"
// #include "mpi.h"
// #include "MeshHandler.hpp"
// #include "CnP.hpp"
// //#include "CnE.hpp"
// //#include "SimulationLoop.hpp"
// #include <iostream>

// using namespace mfem;
// using namespace std;


// int main(int argc, char *argv[])
// {
//     // Initialize MPI and HYPRE.
//     Mpi::Init(argc, argv);
//     Hypre::Init();

//     int myid = Mpi::WorldRank();

//     const char *mesh_file = "Mesh_3x90_T3.mesh";
//     const char *dsF_file = "dsF_3x90_T3.txt";
//     int order = 1;

//     MeshHandler mesh_handler(mesh_file, dsF_file, order);
//     mesh_handler.InitializeMesh();
//     mesh_handler.PrintMeshInfo();

//     InitializeCnP(*mesh_handler.GetFESpace(), *mesh_handler.GetPsi(), mesh_handler.GetGtPsi(), *mesh_handler.GetPmesh());
//     //InitializeCnE(*mesh_handler.GetFESpace(), *mesh_handler.GetPse(), *mesh_handler.GetPmesh());

//     //SimulationLoop loop(mesh_handler);
//     // Inside For Loop
//     // CnP
//     // CnE ...
//     // Reaction
//     //loop.Run(100); 

//     Mpi::Finalize();
//     return 0;
// }

#include "mfem.hpp"
#include "mpi.h"
#include <iostream>
#include "MeshHandler.hpp"
#include "CnP.hpp"
#include "Constants.hpp"

using namespace mfem;
using namespace std;

int main(int argc, char *argv[]) {
    
    // Initialize MPI and HYPRE
    Mpi::Init(argc, argv);
    Hypre::Init();

    // Create the MeshHandler object
    MeshHandler mesh_handler;
    mesh_handler.InitializeMesh();
    mesh_handler.PrintMeshInfo();

    // Create the CnP object
    CnP cnp(mesh_handler);
    cnp.Initialize();

    // Time-stepping loop
    for (int t = 0; t < 10 + 1; ++t) {
        cnp.TimeStep(Constants::dt);
    }

    // Save the results
    cnp.Save();
    
    // Finalize MPI
    Mpi::Finalize();
    return 0;
}
























// // next sections

// // #include "CnE.hpp"
// // #include "PotP.hpp"
// // #include "PotE.hpp"
// // #include "TS_CnP.hpp"

//     // InitializeCnE(*mesh_handler.GetFESpace(), *mesh_handler.GetPse(), *mesh_handler.GetPmesh());
//     // InitializePotP(*mesh_handler.GetFESpace());
//     // InitializePotE(*mesh_handler.GetFESpace());

//     // cout << "BvP: " << BvP << endl;
//     // cout << "BvE: " << BvE << endl;

//     // double Vcell = BvP - BvE;
//     // cout << "Vcell: " << Vcell << endl;