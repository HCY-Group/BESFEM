#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include <fstream>
#include <iostream>


using namespace mfem;
using namespace std;

Concentrations::Concentrations(MeshHandler &mesh_handler)
    : fespace(mesh_handler.GetFESpace()), psi(*mesh_handler.GetPsi()), pse(*mesh_handler.GetPse()),
      EVol(mesh_handler.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi())
{
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();
    
    // Initialize ParGridFunction for CnP and CnE
    CnP = make_unique<ParGridFunction>(fespace);
    CnE = make_unique<ParGridFunction>(fespace);
}


void Concentrations::InitializeCnP() {

    Lithiation(*CnP, 0.3);
    SetupBoundaryConditions();
    SBM_Matrix(psi, Mmatp);
    Solver(Mmatp);

}


// void Concentrations::InitializeCnE() {

//     Lithiation(*CnP, 0.3);


// }


void Concentrations::Lithiation(mfem::ParGridFunction &Cn, double initial_value) {

    Cn = initial_value;

    ParGridFunction TmpF(fespace);
    TmpF = Cn;
    TmpF *= psi;

    double lSum = 0.0;
    Array<double> VtxVal(nC);
    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        TmpF.GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        lSum += EAvg(ei) * EVol(ei);
    }

    double gSum;
    MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double Xfr = gSum / gtPsi;

    // cout << "gSum: " << gSum << std::endl;

}

void Concentrations::SBM_Matrix(mfem::ParGridFunction &psx, HypreParMatrix &Mmat) {

    SetupBoundaryConditions();

    std::unique_ptr<ParBilinearForm> M(new ParBilinearForm(fespace));
    GridFunctionCoefficient cP(&psx);
    M->AddDomainIntegrator(new MassIntegrator(cP));
    M->Assemble();
    
    // HypreParMatrix Mmat;
    M->FormSystemMatrix(boundary_dofs, Mmat); 

    // Mmat.Print("Mmat.txt", 0);


}

void Concentrations::Solver(HypreParMatrix &Mmat) {
    
    HypreSmoother M_prec;
    CGSolver M_solver(MPI_COMM_WORLD);

    M_solver.iterative_mode = false;
    M_solver.SetRelTol(1e-7);
    M_solver.SetAbsTol(0);
    M_solver.SetMaxIter(102);
    M_solver.SetPrintLevel(0);
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(Mmat);

}

void Concentrations::SetupBoundaryConditions() {
    
    // Boundary attributes for Neumann BC on the west boundary
    Array<int> nbc_w_bdr(fespace->GetMesh()->bdr_attributes.Max());
    nbc_w_bdr = 0; 
    nbc_w_bdr[0] = 1;  // Applying Neumann BC to the west boundary

    // Dirichlet BC on the east boundary for CnP
    Array<int> dbc_e_bdr(fespace->GetMesh()->bdr_attributes.Max());
    dbc_e_bdr = 0; 
    dbc_e_bdr[2] = 1;  // Applying Dirichlet BC to the east boundary

    // Extract essential true DOFs (Dirichlet BCs) on the east boundary
    Array<int> ess_tdof_list_e(0);
    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);
    
    // Dirichlet BC on the west boundary for CnE
    Array<int> dbc_w_bdr(fespace->GetMesh()->bdr_attributes.Max());
    dbc_w_bdr = 0; 
    dbc_w_bdr[0] = 1;  // Applying Dirichlet BC to the west boundary

    // Extract essential true DOFs (Dirichlet BCs) on the west boundary
    Array<int> ess_tdof_list_w(0);
    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);
    
    // // Combine the Dirichlet BCs into one array (if needed later)
    // boundary_dofs.SetSize(ess_tdof_list_e.Size() + ess_tdof_list_w.Size());
    // boundary_dofs.CopyFrom(ess_tdof_list_e);
    // boundary_dofs.Append(ess_tdof_list_w);
    
}

