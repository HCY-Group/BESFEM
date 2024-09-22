#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>

// reads mesh_file & dsF_file from Constants.cpp

using namespace mfem;
using namespace std;

MeshHandler::MeshHandler()
    : mesh_file(Constants::mesh_file), dsF_file(Constants::dsF_file), order(Constants::order), dh(Constants::dh), 
      zeta(Constants::zeta), eps(Constants::eps), rho(Constants::rho), Cr(Constants::Cr),
      gtPsi(0.0), gtPse(0.0), gTrgI(0.0) {}

// Function to Initialize and Print Mesh Information
void MeshHandler::LoadMesh() {

    InitializeMesh();
    std::cout << "MeshHandler - fespace initialized with " << fespace->GetNE() << " elements." << std::endl;


    PrintMeshInfo();

    MPI_Barrier(MPI_COMM_WORLD);

}

// Function to Initialize Mesh 
void MeshHandler::InitializeMesh() {
    
    Mesh gmesh(mesh_file);
    gmesh.EnsureNCMesh(true);

    // Create global FE space for distance function.
    H1_FECollection gFec(order, gmesh.Dimension());
    gFespace = make_unique<FiniteElementSpace>(&gmesh, &gFec);

    ReadGlobalDistanceFunction(gFespace); 

    // west boundary size
    Vector Rmin, Rmax;
    gmesh.GetBoundingBox(Rmin, Rmax);
    L_w = Rmax(1) - Rmin(1);

    // local (parallel) mesh
    pmesh = make_unique<ParMesh>(MPI_COMM_WORLD, gmesh);

    nV = pmesh->GetNV();					// number of vertices
	nE = pmesh->GetNE();					// number of elements
	nC = pow(2, pmesh->Dimension());		// number of corner vertices

    Array<double> VtxVal(nC);
    Array<int> gVTX(nC);                    // global indices of corner vertices
    Array<int> VTX(nC);                     // local indices of corner vertices

    // Vector EVol;
    // CalculateElementVolume(nE, pmesh, EVol);
    // // cout << "Size of EVol: " << EVol.Size() << endl; // Debug print to check size

    Vector EVolTemp;
    CalculateElementVolume(nE, pmesh, EVolTemp);
    EVol = EVolTemp;

    // Create Local FE Space
    fespace = make_unique<ParFiniteElementSpace>(pmesh.get(), new H1_FECollection(order, pmesh->Dimension()));
    
    // Local (parallel) GridFunction
    dsF = make_unique<ParGridFunction>(fespace.get());

    // Map local to global element indices
    Array<HYPRE_BigInt> E_L2G;
    pmesh->GetGlobalElementIndices(E_L2G);

    // Map local distance function from global one
    for (int ei = 0; ei < nE; ei++) {
        int gei = E_L2G[ei];
        gmesh.GetElementVertices(gei, gVTX);
        pmesh->GetElementVertices(ei, VTX);
        for (int vi = 0; vi < nC; vi++) {
            (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]);
        }
    }

    InterpolateDomainParameters(nV, fespace);
    CalculateTotalPsi(nV, nE, nC, EVol);
    CalculateTotalPse(nV, nE, nC, EVol);
}

// Calculate Element Volume Function
void MeshHandler::CalculateElementVolume(int nE, const std::unique_ptr<mfem::ParMesh>& pmesh, mfem::Vector& EVol) {
    EVol.SetSize(nE);
	for (int ei = 0; ei < nE; ei++){
		EVol(ei) = pmesh->GetElementVolume(ei);	
	} 
}

// Read global distance function
void MeshHandler::ReadGlobalDistanceFunction(const std::unique_ptr<mfem::FiniteElementSpace>& fespace) {
    gDsF = make_unique<GridFunction>(fespace.get());
    Onm = gDsF->Size();
    ifstream myfile(dsF_file);
    for (int gi = 0; gi < Onm; gi++) {
        myfile >> (*gDsF)(gi);
    }
    myfile.close();
}

void MeshHandler::InterpolateDomainParameters(int nV, const std::unique_ptr<mfem::ParFiniteElementSpace>& fespace) {
    psi = make_unique<ParGridFunction>(fespace.get()); // psi is a pointer here
    pse = make_unique<ParGridFunction>(fespace.get());
    AvP = make_unique<ParGridFunction>(fespace.get());
    AvB = make_unique<ParGridFunction>(fespace.get());

    // interpolate domain parameter from distance function
    for (int vi = 0; vi < nV; vi++) {
        (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (zeta * dh))); // must use * to dereference to access
        (*pse)(vi) = 1.0 - (*psi)(vi);
        (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (zeta * dh)), 2) - 1.0) / (2 * zeta * dh);

        if ((*psi)(vi) < eps) { (*psi)(vi) = eps; }
        if ((*pse)(vi) < eps) { (*pse)(vi) = eps; }
    }

    AvB = std::make_unique<mfem::ParGridFunction>(*AvP);

    for (int vi = 0; vi < nV; vi++) {
        if ((*AvP)(vi) * dh < 1.0e-3) { (*AvP)(vi) = 0.0; }
        if ((*AvB)(vi) * dh < 1.0e-3) { (*AvB)(vi) = 0.0; }
    }
}

// Common Function Used in Pse & Psi
void MeshHandler::CalculateTotals(const std::unique_ptr<mfem::ParGridFunction>& GridFunction, int nV, int nE, int nC, const mfem::Vector& EVol, double& localTotal, double& globalTotal) {
    localTotal = 0.0;

    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Array<double> VtxVal(nC);
        GridFunction->GetNodalValues(ei, VtxVal);
        val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        localTotal += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&localTotal, &globalTotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
}

void MeshHandler::CalculateTotalPsi(int nV, int nE, int nC, const mfem::Vector& EVol) {
    CalculateTotals(psi, nV, nE, nC, EVol, tPsi, gtPsi);
    CalculateTargetCurrent(tPsi);
}

void MeshHandler::CalculateTotalPse(int nV, int nE, int nC, const mfem::Vector& EVol) {
    CalculateTotals(pse, nV, nE, nC, EVol, tPse, gtPse);
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

const mfem::Vector& MeshHandler::GetElementVolume() const {
    return EVol; // Ensure this is a reference to the actual Vector
}

void MeshHandler::TestFESpace() {
    // Create a test ParGridFunction on the fespace
    ParGridFunction test_f(fespace.get());
    test_f = 1.0; // Set all values in test_f to 1

    // Output the first value of test_f as a test
    std::cout << "MeshHandler - First value of test_f: " << test_f(0) << std::endl;
}
