#include "MeshHandler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

MeshHandler::MeshHandler()
    : mesh_file(Constants::mesh_file), dsF_file(Constants::dsF_file), order(Constants::order), dh(Constants::dh), 
      zeta(Constants::zeta), eps(Constants::eps), rho(Constants::rho), Cr(Constants::Cr),
      gmesh(mesh_file), gtPsi(0.0), gtPse(0.0), gTrgI(0.0) {}

void MeshHandler::LoadMesh() {
    InitializeMesh();
    PrintMeshInfo();
}

void MeshHandler::InitializeMesh() {
    gmesh.EnsureNCMesh(true);
    gFespace = make_unique<FiniteElementSpace>(&gmesh, new H1_FECollection(order, gmesh.Dimension()));

    ReadGlobalDistanceFunction();

    pmesh = make_unique<ParMesh>(MPI_COMM_WORLD, gmesh);
    fespace = make_unique<ParFiniteElementSpace>(pmesh.get(), new H1_FECollection(order, pmesh->Dimension()));
    dsF = make_unique<ParGridFunction>(fespace.get());

    // Map local distance function from global one
    Array<HYPRE_BigInt> E_L2G;
    pmesh->GetGlobalElementIndices(E_L2G);

    int nV = pmesh->GetNV();					// number of vertices
	int nE = pmesh->GetNE();					// number of elements
	int nC = pow(2, pmesh->Dimension());		// number of corner vertices

    Array<int> gVTX(nC);
    Array<int> VTX(nC);

    for (int ei = 0; ei < nE; ei++) {
        int gei = E_L2G[ei];
        gmesh.GetElementVertices(gei, gVTX);
        pmesh->GetElementVertices(ei, VTX);
        for (int vi = 0; vi < nC; vi++) {
            (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]);
        }
    }

    InterpolateDomainParameters();
    CalculateTotalPsi();
    CalculateTotalPse();
}

void MeshHandler::ReadGlobalDistanceFunction() {
    gDsF = make_unique<GridFunction>(gFespace.get());
    ifstream myfile(dsF_file);
    for (int gi = 0; gi < gDsF->Size(); gi++) {
        myfile >> (*gDsF)(gi);
    }
    myfile.close();
}

void MeshHandler::InterpolateDomainParameters() {
    psi = make_unique<ParGridFunction>(fespace.get());
    pse = make_unique<ParGridFunction>(fespace.get());
    AvP = make_unique<ParGridFunction>(fespace.get());

    for (int vi = 0; vi < fespace->GetNV(); vi++) {
        (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (zeta * dh)));
        (*pse)(vi) = 1.0 - (*psi)(vi);
        (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (zeta * dh)), 2) - 1.0) / (2 * zeta * dh);

        if ((*psi)(vi) < eps) { (*psi)(vi) = eps; }
        if ((*pse)(vi) < eps) { (*pse)(vi) = eps; }
    }

    for (int vi = 0; vi < fespace->GetNV(); vi++) {
        if ((*AvP)(vi) * dh < 1.0e-3) { (*AvP)(vi) = 0.0; }
    }
}

void MeshHandler::CalculateTotalPsi() {
    double tPsi = 0.0;
    
    int nV = pmesh->GetNV();					// number of vertices
	int nE = pmesh->GetNE();					// number of elements
	int nC = pow(2, pmesh->Dimension());		// number of corner vertices

    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Array<double> VtxVal(nC);
        psi->GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        tPsi += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&tPsi, &gtPsi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    CalculateTargetCurrent(tPsi);
}

void MeshHandler::CalculateTotalPse() {
    double tPse = 0.0;

    int nV = pmesh->GetNV();					// number of vertices
	int nE = pmesh->GetNE();					// number of elements
	int nC = pow(2, pmesh->Dimension());		// number of corner vertices

    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Array<double> VtxVal(nC);
        pse->GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        tPse += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&tPse, &gtPse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MeshHandler::CalculateTargetCurrent(double tPsi) {
    double trgI = tPsi * rho * (0.9 - 0.3) / (3600.0 / Cr);
    MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MeshHandler::PrintMeshInfo() {
    
    int nV = pmesh->GetNV();					// number of vertices
	int nE = pmesh->GetNE();					// number of elements
	int nC = pow(2, pmesh->Dimension());		// number of corner vertices
    
    cout << "Number of vertices: " << nV << endl;
    cout << "Number of elements: " << nE << endl;

    Vector Rmin, Rmax;
    gmesh.GetBoundingBox(Rmin, Rmax);
    double L_w = Rmax(1) - Rmin(1);
    cout << "West boundary size: " << L_w << endl;

    cout << "Total Psi: " << gtPsi << endl;
    cout << "Total Pse: " << gtPse << endl;
    cout << "Target Current: " << gTrgI << endl;

    // // psi info

    // cout << "PSI in MeshHandler" << endl;
    // psi->Print(std::cout);
    // cout << "end" << endl;

}

void MeshHandler::Save() {
    if (pmesh) {
        pmesh->Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/Pmesh");
    } else {
        std::cerr << "Error: pmesh is not initialized.\n";
    }
}