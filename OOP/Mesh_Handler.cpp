#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>

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

    // Create global FE space for distance function.
    mfem::H1_FECollection gFec(order, gmesh.Dimension());
    gFespace = make_unique<mfem::FiniteElementSpace>(&gmesh, &gFec);

    ReadGlobalDistanceFunction(gFespace); 

    // west boundary size
    mfem::Vector Rmin, Rmax;
    gmesh.GetBoundingBox(Rmin, Rmax);
    L_w = Rmax(1) - Rmin(1);

    // local (parallel) mesh
    pmesh0 = make_unique<mfem::ParMesh>(MPI_COMM_WORLD, gmesh);

    nV = pmesh0->GetNV();					// number of vertices
	nE = pmesh0->GetNE();					// number of elements
	nC = pow(2, pmesh0->Dimension());		// number of corner vertices

    mfem::Array<double> VtxVal(nC);               // values of corner vertices
    mfem::Array<int> gVTX(nC);                    // global indices of corner vertices
    mfem::Array<int> VTX(nC);                     // local indices of corner vertices

    mfem::Vector EVolTemp;
    CalculateElementVolume(pmesh0, EVolTemp);
    EVol = EVolTemp;

    // Create Local FE Space
    fespace = make_shared<mfem::ParFiniteElementSpace>(pmesh0.get(), new mfem::H1_FECollection(order, pmesh0->Dimension()));

    // Local (parallel) GridFunction
    dsF = make_unique<mfem::ParGridFunction>(fespace.get());

    // Map local to global element indices
    mfem::Array<HYPRE_BigInt> E_L2G;
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

    InitializeGridFunctions(fespace);
    InterpolateDomainParameters(fespace);
    CalculatePhasePotentialsAndTargetCurrent();
}

// Calculate Element Volume Function
void MeshHandler::CalculateElementVolume(const std::unique_ptr<mfem::ParMesh>& pmesh, mfem::Vector& EVol) {
    EVol.SetSize(nE);
	for (int ei = 0; ei < nE; ei++){
        EVol(ei) = pmesh->GetElementVolume(ei);	
	} 
}

// Read global distance function
void MeshHandler::ReadGlobalDistanceFunction(const std::unique_ptr<mfem::FiniteElementSpace>& fespace) {
    gDsF = make_unique<mfem::GridFunction>(fespace.get());
    Onm = gDsF->Size();
    ifstream myfile(dsF_file);
    for (int gi = 0; gi < Onm; gi++) {
        myfile >> (*gDsF)(gi);
    }
    myfile.close();
}

void MeshHandler::InitializeGridFunctions(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace) {
    psi = make_unique<mfem::ParGridFunction>(fespace.get());
    pse = make_unique<mfem::ParGridFunction>(fespace.get());
    AvP = make_unique<mfem::ParGridFunction>(fespace.get());
    AvB = make_unique<mfem::ParGridFunction>(fespace.get());
}

void MeshHandler::InterpolateDomainParameters(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace) {

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

void MeshHandler::CalculateTotals(const mfem::ParGridFunction& grid_function, const mfem::Vector& element_volumes, double& local_total, double& global_total) {
    local_total = 0.0;

    // Average value for each element
    mfem::Vector element_avg_values(nE);

    for (int ei = 0; ei < nE; ei++) {
        mfem::Array<double> vertex_values(nC);
        grid_function.GetNodalValues(ei, vertex_values);

        // Compute the average value over all corners of the element
        double avg_value = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            avg_value += vertex_values[vt];
        }
        avg_value /= nC;

        element_avg_values(ei) = avg_value;

        // Accumulate the weighted value for the current element
        local_total += avg_value * element_volumes(ei);
    }

    // Perform global MPI reduction to sum values across all processes
    MPI_Allreduce(&local_total, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MeshHandler::CalculateTotalPhaseField(const mfem::ParGridFunction& grid_function, double& total, double& global_total) {
    
    // Call the general `CalculateTotals` method for any field
    CalculateTotals(grid_function, EVol, total, global_total);

}

void MeshHandler::CalculatePhasePotentialsAndTargetCurrent() {
    
    // Calculate totals for Psi and Pse fields
    CalculateTotalPhaseField(*psi, tPsi, gtPsi);
    CalculateTotalPhaseField(*pse, tPse, gtPse);

    // Compute the target current using Psi
    CalculateTargetCurrent(tPsi);
}

void MeshHandler::CalculateTargetCurrent(double total_psi) {
    
    // Compute target current based on total Psi, rho, Cr, and constants
    trgI = total_psi * rho * (0.9 - 0.3) / (3600.0 / Cr);

    // Perform global MPI reduction to get the total target current
    MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void MeshHandler::PrintMeshInfo() {
    
    cout << "Number of vertices: " << nV << endl;
    cout << "Number of elements: " << nE << endl;

    cout << "West boundary size: " << L_w << endl;

    cout << "Total Psi: " << gtPsi << endl;
    cout << "Total Pse: " << gtPse << endl;
    cout << "Target Current: " << gTrgI << endl;

}

void MeshHandler::SetupBoundaryConditions(mfem::ParMesh *pmesh, mfem::ParFiniteElementSpace *fespace) {
        
    // Boundary attributes for Neumann BC on the west boundary
    nbc_w_bdr.SetSize(pmesh->bdr_attributes.Max());
    nbc_w_bdr = 0; 
    nbc_w_bdr[0] = 1;  // Applying Neumann BC to the west boundary

    // Dirichlet BC on the east boundary for CnP
    dbc_e_bdr.SetSize(pmesh->bdr_attributes.Max());
    dbc_e_bdr = 0; 
    dbc_e_bdr[2] = 1;  // Applying Dirichlet BC to the east boundary

    // Extract essential true DOFs (Dirichlet BCs) on the east boundary
    // mfem::Array<int> ess_tdof_list_e(0);
    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

    // Dirichlet BC on the west boundary for CnE
    dbc_w_bdr.SetSize(pmesh->bdr_attributes.Max());
    dbc_w_bdr = 0; 
    dbc_w_bdr[0] = 1;  // Applying Dirichlet BC to the west boundary

    // Extract essential true DOFs (Dirichlet BCs) on the west boundary
    // mfem::Array<int> ess_tdof_list_w(0);
    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);
    
}

