// mpicxx -g -O3 -std=c++14 -I../../.. -I../../../../hypre/src/hypre/include main2.cpp 
// Constants.cpp MeshHandler.cpp CnP.cpp CnE.cpp PotP.cpp PotE.cpp Reaction.cpp -o main2 
// -L../../.. -lmfem -L../../../../hypre/src/hypre/lib
// -lHYPRE -L../../../../metis-4.0 -lmetis -lrt

#include "mfem.hpp"
#include "mpi.h"
#include <iostream>
#include "MeshHandler.hpp"
#include "CnP.hpp"
#include "CnE.hpp"
#include "Constants.hpp"
#include "PotP.hpp"
#include "PotE.hpp"
#include "Reaction.hpp"

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
    //mesh_handler.Save();

    // Create the CnP object
    CnP cnp(mesh_handler);
    cnp.Initialize();

    // Create the CnE object
    CnE cne(mesh_handler, cnp);
    cne.Initialize();

    // Create the PotP object
    PotP potp(mesh_handler);
    potp.Initialize();

    // Create the PotE object
    PotE pote(mesh_handler, potp);
    pote.Initialize();

    // Create Reaction
    // Reaction reaction(mesh_handler, cne);


    //Time-stepping loop
    for (int t = 0; t < 10 + 1; ++t) {
        cnp.TimeStep(Constants::dt);
        cne.TimeStep(Constants::dt);
        
        // Create Reaction
        Reaction reaction(mesh_handler, cne, pote);
        reaction.TimeStep(Constants::dt);


    }




    // // Save the results
    // cnp.Save();
    // cne.Save();


    // Finalize MPI
    Mpi::Finalize();
    return 0;
}
