#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include <fstream>
#include <iostream>


using namespace mfem;
using namespace std;

Concentrations::Concentrations(MeshHandler &mesh_handler)
    : fespace(mesh_handler.GetFESpace()), psi(*mesh_handler.GetPsi()), pse(*mesh_handler.GetPse()),
      EVol(mesh_handler.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi()), reaction(mesh_handler, *this)
      
{
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();
    
    // Initialize ParGridFunction for CnP and CnE
    CnP = make_unique<ParGridFunction>(fespace);
    CnE = make_unique<ParGridFunction>(fespace);

    Rxc = make_unique<ParGridFunction>(fespace);
    Rxe = make_unique<ParGridFunction>(fespace);


    reaction.Initialize();
}


void Concentrations::InitializeCnP() {

    Lithiation(*CnP, 0.3);
    SetupBoundaryConditions();
    SBM_Matrix(psi, Mmatp);
    CGSolver Mp_solver(MPI_COMM_WORLD);
    Solver(Mmatp, Mp_solver); 

}

void Concentrations::InitializeCnE() {

    CreateCnE(*CnE, 0.001);
    SetupBoundaryConditions();
    SBM_Matrix(pse, Mmate);
    CGSolver Me_solver(MPI_COMM_WORLD);
    Solver(Mmate, Me_solver);

}

void Concentrations::TimeStepCnP() {

    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    GridFunctionCoefficient cAp(Rxc.get());
    SetupRx(Rxn, *Rxc, Constants::rho, cAp); // Rxn needs to come from Reacions.cpp

    // Rxc->Print(std::cout);

    // cout << "rho : " << Constants::rho << endl;

}

void Concentrations::TimeStepCnE() {

    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    GridFunctionCoefficient cAe(Rxe.get());
    SetupRx(Rxn, *Rxe, Constants::t_minus, cAe);

    Rxe->Print(std::cout);

}

void Concentrations::CreateCnE(mfem::ParGridFunction &Cn, double initial_value) {

    Cn = initial_value;

}

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
    
    M->FormSystemMatrix(boundary_dofs, Mmat); 

}

void Concentrations::Solver(HypreParMatrix &Mmat, CGSolver &M_solver) {
    
    HypreSmoother M_prec;

    M_solver.iterative_mode = false;
    M_solver.SetRelTol(1e-7);
    M_solver.SetAbsTol(0);
    M_solver.SetMaxIter(102);
    M_solver.SetPrintLevel(0);
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.SetPreconditioner(M_prec);

    if (&Mmat == &Mmatp) {
        M_solver.SetOperator(Mmat); // this is only needed for Mmatp
    }

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
    
}

void Concentrations::SetupRx(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, 
double value, GridFunctionCoefficient cAx) {

    Rx2 = Rx1;

    if (&Rx2 == Rxc.get()) {
        Rx2 /= value; // this is needed for CnP & Rxc
    }

    if (&Rx2 == Rxe.get()) {
        Rx2 *= (-1.0 * value); // this is needed for CnE & Rxe
    }

}
