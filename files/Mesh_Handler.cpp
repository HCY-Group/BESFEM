#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>

// reads mesh_file & dsF_file from Constants.cpp

using namespace mfem;
using namespace std;

double gTrgI = 0.0;


MeshHandler::MeshHandler()
    : mesh_file(Constants::mesh_file), dsF_file(Constants::dsF_file), order(Constants::order), dh(Constants::dh), 
      zeta(Constants::zeta), eps(Constants::eps), rho(Constants::rho), Cr(Constants::Cr),
      gtPsi(0.0), gtPse(0.0){}


// Function to Initialize and Print Mesh Information
void MeshHandler::LoadMesh() {

    InitializeMesh();
    PrintMeshInfo();

    MPI_Barrier(MPI_COMM_WORLD);

}

// Function to Initialize Mesh 
void MeshHandler::InitializeMesh() {
    
    Mesh gmesh(mesh_file);
    gmesh.EnsureNCMesh(true);

    // // Uniformly refine the mesh, e.g., 2 times
    // int number_of_refinements = 2; // Change as needed
    // for (int i = 0; i < number_of_refinements; i++) {
    //     gmesh.UniformRefinement();
    // }

    // Create global FE space for distance function.
    H1_FECollection gFec(order, gmesh.Dimension());
    gFespace = make_unique<FiniteElementSpace>(&gmesh, &gFec);

    ReadGlobalDistanceFunction(gFespace); 

    // west boundary size
    Vector Rmin, Rmax;
    gmesh.GetBoundingBox(Rmin, Rmax);
    L_w = Rmax(1) - Rmin(1);

    // local (parallel) mesh
    pmesh0 = make_unique<ParMesh>(MPI_COMM_WORLD, gmesh);

    nV = pmesh0->GetNV();					// number of vertices
	nE = pmesh0->GetNE();					// number of elements
	nC = pow(2, pmesh0->Dimension());		// number of corner vertices

    Array<double> VtxVal(nC);
    Array<int> gVTX(nC);                    // global indices of corner vertices
    Array<int> VTX(nC);                     // local indices of corner vertices

    // Vector EVol;
    // CalculateElementVolume(nE, pmesh, EVol);
    // // cout << "Size of EVol: " << EVol.Size() << endl; // Debug print to check size

    Vector EVolTemp;
    CalculateElementVolume(nE, pmesh0, EVolTemp);
    EVol = EVolTemp;

    // Create Local FE Space
    fespace = make_shared<ParFiniteElementSpace>(pmesh0.get(), new H1_FECollection(order, pmesh0->Dimension()));
    // H1_FECollection fec(order, pmesh->Dimension());	
    // fespace = std::make_shared<mfem::ParFiniteElementSpace>(pmesh, fec);

    

    // Local (parallel) GridFunction
    dsF = make_unique<ParGridFunction>(fespace.get());

    // Map local to global element indices
    Array<HYPRE_BigInt> E_L2G;
    pmesh0->GetGlobalElementIndices(E_L2G);

    // Map local distance function from global one
    for (int ei = 0; ei < nE; ei++) {
        int gei = E_L2G[ei];
        gmesh.GetElementVertices(gei, gVTX);
        pmesh0->GetElementVertices(ei, VTX);
        for (int vi = 0; vi < nC; vi++) {
            (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]);
        }
    }

    pmesh = std::move(pmesh0);

    InterpolateDomainParameters(nV, fespace);
    CalculateTotalPsi(nV, nE, nC, EVol);
    CalculateTotalPse(nV, nE, nC, EVol);
    // SetupBoundaryConditions(pmesh.get(), fespace.get());
}

mfem::ParMesh MeshHandler::GetMesh() {
    mfem::ParMesh tmpmesh(*pmesh);
    return std::move(tmpmesh);
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

void MeshHandler::InterpolateDomainParameters(int nV, const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace) {
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
        if ((*AvB)(vi) * dh < 1.0e-6) { (*AvB)(vi) = 0.0; }
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

void MeshHandler::SetupBoundaryConditions(mfem::ParMesh *pmesh, mfem::ParFiniteElementSpace *fespace) {
    
    Array<int> boundary_dofs;
    
    // Boundary attributes for Neumann BC on the west boundary
    Array<int> nbc_w_bdr(pmesh->bdr_attributes.Max());
    nbc_w_bdr = 0; 
    nbc_w_bdr[0] = 1;  // Applying Neumann BC to the west boundary

    // // Printing the values
    // std::cout << "nbc_w_bdr values: ";
    // for (int i = 0; i < nbc_w_bdr.Size(); i++) {
    //     std::cout << nbc_w_bdr[i] << " ";
    // }
    // std::cout << std::endl;

    // Dirichlet BC on the east boundary for CnP
    Array<int> dbc_e_bdr(pmesh->bdr_attributes.Max());
    dbc_e_bdr = 0; 
    dbc_e_bdr[2] = 1;  // Applying Dirichlet BC to the east boundary

    // std::cout << "dbc_e_bdr size in MeshHandler: " << dbc_e_bdr.Size() << std::endl;


    // Extract essential true DOFs (Dirichlet BCs) on the east boundary
    Array<int> ess_tdof_list_e(0);
    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

    // Dirichlet BC on the west boundary for CnE
    Array<int> dbc_w_bdr(pmesh->bdr_attributes.Max());
    dbc_w_bdr = 0; 
    dbc_w_bdr[0] = 1;  // Applying Dirichlet BC to the west boundary

    // Extract essential true DOFs (Dirichlet BCs) on the west boundary
    Array<int> ess_tdof_list_w(0);
    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);
    
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

std::shared_ptr<mfem::ParFiniteElementSpace> MeshHandler::GetFESpace() {
    return fespace;  // Provide access to the fespace
}

