#include "mfem.hpp"
#include "mpi.h"
#include "MeshHandler.hpp"
#include "CnP.hpp"
#include "CnE.hpp"
#include "PotP.hpp"
#include "PotE.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
    // Initialize MPI and HYPRE.
    Mpi::Init(argc, argv);
    Hypre::Init();

    int myid = Mpi::WorldRank();

    const char *mesh_file = "Mesh_3x90_T3.mesh";
    const char *dsF_file = "dsF_3x90_T3.txt";
    int order = 1;

    MeshHandler mesh_handler(mesh_file, dsF_file, order);
    mesh_handler.InitializeMesh();
    mesh_handler.PrintMeshInfo();

    InitializeCnP(*mesh_handler.GetFESpace(), *mesh_handler.GetPsi(), mesh_handler.GetGtPsi(), *mesh_handler.GetPmesh());
    InitializeCnE(*mesh_handler.GetFESpace(), *mesh_handler.GetPse(), *mesh_handler.GetPmesh());
    InitializePotP(*mesh_handler.GetFESpace());
    InitializePotE(*mesh_handler.GetFESpace());

    // cout << "BvP: " << BvP << endl;
    // cout << "BvE: " << BvE << endl;

    double Vcell = BvP - BvE;
    cout << "Vcell: " << Vcell << endl;


    Mpi::Finalize();
    return 0;
}


