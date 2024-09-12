#include "MeshHandler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>

// reads mesh_file & dsF_file from Constants.cpp

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

    // Create global FE space for distance function.
    gFespace = make_unique<FiniteElementSpace>(&gmesh, new H1_FECollection(order, gmesh.Dimension()));

    ReadGlobalDistanceFunction();

    // west boundary size
    Vector Rmin, Rmax;
    gmesh.GetBoundingBox(Rmin, Rmax);
    L_w = Rmax(1) - Rmin(1);

    // local (parallel) mesh
    pmesh = make_unique<ParMesh>(MPI_COMM_WORLD, gmesh);
    fespace = make_unique<ParFiniteElementSpace>(pmesh.get(), new H1_FECollection(order, pmesh->Dimension()));
    dsF = make_unique<ParGridFunction>(fespace.get());

    // Map local distance function from global one
    Array<HYPRE_BigInt> E_L2G;
    pmesh->GetGlobalElementIndices(E_L2G);

    nV = pmesh->GetNV();					// number of vertices
	nE = pmesh->GetNE();					// number of elements
	nC = pow(2, pmesh->Dimension());		// number of corner vertices

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

// Read global distance function
void MeshHandler::ReadGlobalDistanceFunction() {
    gDsF = make_unique<GridFunction>(gFespace.get());
    ifstream myfile(dsF_file);
    for (int gi = 0; gi < gDsF->Size(); gi++) {
        myfile >> (*gDsF)(gi);
    }
    myfile.close();
}

void MeshHandler::InterpolateDomainParameters() {
    psi = make_unique<ParGridFunction>(fespace.get()); // psi is a pointer here
    pse = make_unique<ParGridFunction>(fespace.get());
    AvP = make_unique<ParGridFunction>(fespace.get());
    AvB = make_unique<ParGridFunction>(fespace.get());

    for (int vi = 0; vi < fespace->GetNV(); vi++) {
        (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (zeta * dh))); // must use * to dereference to access
        (*pse)(vi) = 1.0 - (*psi)(vi);
        (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (zeta * dh)), 2) - 1.0) / (2 * zeta * dh);

        if ((*psi)(vi) < eps) { (*psi)(vi) = eps; }
        if ((*pse)(vi) < eps) { (*pse)(vi) = eps; }
    }

    AvB = std::make_unique<mfem::ParGridFunction>(*AvP);

    for (int vi = 0; vi < fespace->GetNV(); vi++) {
        if ((*AvP)(vi) * dh < 1.0e-3) { (*AvP)(vi) = 0.0; }
        if ((*AvB)(vi) * dh < 1.0e-3) { (*AvB)(vi) = 0.0; }
    }
}

void MeshHandler::CalculateTotalPsi() {
    tPsi = 0.0;
    
    nV = pmesh->GetNV();					// number of vertices
	nE = pmesh->GetNE();					// number of elements
	nC = pow(2, pmesh->Dimension());		// number of corner vertices

    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Array<double> VtxVal(nC);
        psi->GetNodalValues(ei, VtxVal);
        val = 0.0;
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
    tPse = 0.0;

    nV = pmesh->GetNV();					// number of vertices, from MFEM parmesh
	nE = pmesh->GetNE();					// number of elements, from MFEM parmesh
	nC = pow(2, pmesh->Dimension());		// number of corner vertices

    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Array<double> VtxVal(nC);
        pse->GetNodalValues(ei, VtxVal);
        val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        tPse += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&tPse, &gtPse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MeshHandler::CalculateTargetCurrent(double tPsi) {
    trgI = tPsi * rho * (0.9 - 0.3) / (3600.0 / Cr);
    MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MeshHandler::PrintMeshInfo() {
    
    nV = pmesh->GetNV();					// number of vertices
	nE = pmesh->GetNE();					// number of elements
	nC = pow(2, pmesh->Dimension());		// number of corner vertices
    
    cout << "Number of vertices: " << nV << endl;
    cout << "Number of elements: " << nE << endl;

    cout << "West boundary size: " << L_w << endl;

    cout << "Total Psi: " << gtPsi << endl;
    cout << "Total Pse: " << gtPse << endl;
    cout << "Target Current: " << gTrgI << endl;

}

void MeshHandler::Save() {
    if (pmesh) {
        pmesh->Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/Pmesh");
    } else {
        std::cerr << "Error: pmesh is not initialized.\n";
    }
}