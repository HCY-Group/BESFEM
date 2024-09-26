// mpicxx -g -O3 -std=c++14 -I../../.. -I../../../../hypre/src/hypre/include main2.cpp 
// Constants.cpp MeshHandler.cpp CnP.cpp CnE.cpp PotP.cpp PotE.cpp Reaction.cpp -o main 
// -L../../.. -lmfem -L../../../../hypre/src/hypre/lib
// -lHYPRE -L../../../../metis-4.0 -lmetis -lrt

// mpicxx -g -O3 -std=c++14 -I../../.. -I../../../../hypre/src/hypre/include main.cpp Constants.cpp Mesh_Handler.cpp -o main -L../../.. -lmfem -L../../../../hypre/src/hypre/lib -lHYPRE -L../../../../metis-4.0 -lmetis -lrt

#include "mfem.hpp"
#include "mpi.h"

#include "Constants.hpp"
#include "Mesh_Handler.hpp"
#include "Concentrations.hpp"
#include "Reaction.hpp"

#include <iostream>

using namespace mfem;
using namespace std;

int main(int argc, char *argv[]) {
    
    // Initialize MPI and HYPRE
    Mpi::Init(argc, argv);
    Hypre::Init();

    {

    // Create the MeshHandler object
    MeshHandler mesh_handler;       // define mesh and dsf file in Constants.cpp
    mesh_handler.LoadMesh();
    //mesh_handler.Save();

    // mesh_handler.TestFESpace();

    // Initialize CnP & CnE
    Concentrations concentrations(mesh_handler);
    // concentrations.TestFESpace(mesh_handler.GetFESpace());

    // MFEM_VERIFY(MeshHandler.GetFESpace()->GetMesh() == Concentrations.GetFESpace()->GetMesh(), "Mesh mismatch between modules!");


    
    concentrations.SetupBoundaryConditions(mesh_handler.GetFESpace());

    concentrations.InitializeCnP(mesh_handler.GetFESpace());
    concentrations.InitializeCnE(mesh_handler.GetFESpace());  

    // // Initialize Reaction
    Reaction reaction(mesh_handler, concentrations);
    reaction.Initialize(); 

    //Time-stepping loop
    for (int t = 0; t < 10 + 1; ++t) {
        concentrations.TimeStepCnP(mesh_handler.GetFESpace());
        concentrations.TimeStepCnE(mesh_handler.GetFESpace());


    }

    }




    // // Save the results
    // cnp.Save();
    // cne.Save();


    // Finalize MPI
    Mpi::Finalize();
    return 0;

}