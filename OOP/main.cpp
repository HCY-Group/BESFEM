// mpicxx -g -O3 -std=c++14 -I../../.. -I../../../../hypre/src/hypre/include main2.cpp 
// Constants.cpp MeshHandler.cpp CnP.cpp CnE.cpp PotP.cpp PotE.cpp Reaction.cpp -o main 
// -L../../.. -lmfem -L../../../../hypre/src/hypre/lib
// -lHYPRE -L../../../../metis-4.0 -lmetis -lrt

// mpicxx -g -O3 -std=c++14 -I../../.. -I../../../../hypre/src/hypre/include main.cpp Constants.cpp Mesh_Handler.cpp -o main -L../../.. -lmfem -L../../../../hypre/src/hypre/lib -lHYPRE -L../../../../metis-4.0 -lmetis -lrt

#include "mfem.hpp"
#include "mpi.h"

#include "Constants.hpp"
#include "Mesh_Handler.hpp"

// #include "CnP.hpp"
// #include "CnE.hpp"
// #include "PotP.hpp"
// #include "PotE.hpp"
// #include "Reaction.hpp"

#include <iostream>

using namespace mfem;
using namespace std;

int main(int argc, char *argv[]) {
    
    // Initialize MPI and HYPRE
    Mpi::Init(argc, argv);
    Hypre::Init();

    // Create the MeshHandler object
    MeshHandler mesh_handler;       // define mesh and dsf file in Constants.cpp
    mesh_handler.LoadMesh();
    //mesh_handler.Save();

    // // Create Reaction
    // Reaction reaction(mesh_handler);
    // reaction.Initialize(); 

    // // Create the CnP object
    // CnP cnp(mesh_handler, reaction);
    // cnp.Initialize();

    // // Create the CnE object
    // CnE cne(mesh_handler, reaction, cnp);
    // cne.Initialize();

    // // Create the PotP object
    // PotP potp(mesh_handler);
    // potp.Initialize();

    // // Create the PotE object
    // PotE pote(mesh_handler, potp, cne);
    // pote.Initialize();

    // // // Create Reaction
    // // // Reaction reaction(mesh_handler, cne, pote);

    // //Time-stepping loop
    // for (int t = 0; t < 10 + 1; ++t) {
    //     cnp.TimeStep(Constants::dt);
    //     cne.TimeStep(Constants::dt);
    //     // potp.TimeStep(Constants::dt);
    //     pote.TimeStep(Constants::dt);
    //     // reaction.TimeStep(Constants::dt);
        
    //     // Create Reaction
    //     // Reaction reaction(mesh_handler, cne, pote);
    //     // reaction.TimeStep(Constants::dt);
    //     // for ......
    //     // end


    // }




    // // // Save the results
    // // cnp.Save();
    // // cne.Save();


    // Finalize MPI
    Mpi::Finalize();
    return 0;
}