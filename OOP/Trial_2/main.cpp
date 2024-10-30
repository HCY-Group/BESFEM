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
#include "Potentials.hpp"

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
    // mesh_handler.Save();

    // Create pmesh & fespace to use for all 
    ParMesh pmesh = mesh_handler.GetMesh();
    H1_FECollection fec(Constants::order, pmesh.Dimension());	
	ParFiniteElementSpace fespace(&pmesh, &fec);
    mesh_handler.SetupBoundaryConditions(&pmesh, &fespace); 

    // Initialize CnP & CnE
    Concentrations concentrations(&pmesh, &fespace, mesh_handler);
    concentrations.InitializeCnP();
    concentrations.InitializeCnE();  

    // Initialize Reaction
    // Reaction reaction(&fespace, mesh_handler, concentrations);
    // reaction.Initialize();
    concentrations.reaction.Initialize(); 

    // Initialize PotP & PotE
    // Potentials potentials(&pmesh, &fespace, mesh_handler);
    concentrations.potentials.InitializePotP();
    concentrations.potentials.InitializePotE();

    // Time-stepping loop
    for (int t = 0; t < 10 + 1; ++t) {
        concentrations.TimeStepCnP();
        concentrations.TimeStepCnE();
        concentrations.reaction.TimeStep();

    }

    // Save the Results
    // concentrations.SaveCnP();
    // concentrations.SaveCnE();




    }






    // Finalize MPI
    Mpi::Finalize();
    return 0;

}