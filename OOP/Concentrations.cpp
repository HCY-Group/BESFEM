#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include <fstream>
#include <iostream>


using namespace mfem;
using namespace std;

Concentrations::Concentrations(MeshHandler &mesh_handler)
    : mesh_handler(mesh_handler), fespace(mesh_handler.GetFESpace()), pmesh(mesh_handler.GetPmesh()), psi(*mesh_handler.GetPsi()), pse(*mesh_handler.GetPse()),
      EVol(mesh_handler.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi()), reaction(mesh_handler, *this), PeR(fespace), matCoef_R(&PeR)
      
{

    nbc_w_bdr.SetSize(mesh_handler.pmesh->bdr_attributes.Max());
    nbc_w_bdr = 0;  // Initialize the array with zeros

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

// Concentrations::Concentrations(MeshHandler &mesh_handler)
//     : mesh_handler(mesh_handler), fespace(mesh_handler.GetFESpace()) {
//     std::cout << "Concentrations - fespace pointer: " << fespace.get() << std::endl;

//     mfem::ParGridFunction& myTestF = mesh_handler.GetTestF();
//     std::cout << "First value of test_f in Concentrations: " << myTestF(0) << std::endl; 


//     // ParGridFunction test_f(fespace.get());
//     // std::cout << "Concentrations - Initial value of test_f: " << test_f(0) << std::endl;


//     // Initialize ParGridFunction for CnP and CnE using the shared fespace
//     // CnP = std::make_unique<ParGridFunction>(fespace.get());
//     // CnE = std::make_unique<ParGridFunction>(fespace.get());

//     // Rxc = std::make_unique<ParGridFunction>(fespace.get());
//     // Rxe = std::make_unique<ParGridFunction>(fespace.get());

// }


void Concentrations::InitializeCnP(ParFiniteElementSpace *fespace) {

    Lithiation(*CnP, 0.3, fespace);
    // SetupBoundaryConditions();
    SBM_Matrix(psi, Mmatp, fespace);
    CGSolver Mp_solver(MPI_COMM_WORLD);
    Solver(Mmatp, Mp_solver); 

    // CnP->Print(std::cout);

}

void Concentrations::InitializeCnE(ParFiniteElementSpace *fespace) {

    CreateCnE(*CnE, 0.001);
    // SetupBoundaryConditions();
    SBM_Matrix(pse, Mmate, fespace);
    CGSolver Me_solver(MPI_COMM_WORLD);
    Solver(Mmate, Me_solver);
    ImposeNeumannBC(PeR, pse);
    // GridFunctionCoefficient matCoef_R(&PeR);

    // PeR.Print(std::cout);

}

void Concentrations::TimeStepCnP(ParFiniteElementSpace *fespace) {

    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    GridFunctionCoefficient cAp(Rxc.get());
    SetupRx(Rxn, *Rxc, Constants::rho, cAp); // Rxn needs to come from Reacions.cpp

    // Dummy boundary: empty array
    Array<int> dummy_boundary;
    
    // Dummy coefficient: 0.0 (won't affect the result if not used)
    ConstantCoefficient dummy_coef(0.0);

    // Call ForceTerm with dummy boundary and coefficient
    ForceTerm(fespace, cAp, Fct, dummy_boundary, dummy_coef, false);
    // Fct.Print(std::cout);

    // Rxc->Print(std::cout);



}

void Concentrations::TimeStepCnE(ParFiniteElementSpace *fespace) {

    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    GridFunctionCoefficient cAe(Rxe.get());
    SetupRx(Rxn, *Rxe, Constants::t_minus, cAe);

    TotalReaction(*Rxe, eCrnt);
    ConstantCoefficient nbcCoef(-infx); // Neumann BC works with -infx

    // PeR.Print(std::cout);

    ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    ForceTerm(fespace, cAe, Fet, nbc_w_bdr, nbcCoef, true); // true since applying boundary conditions

    Fet.Print(std::cout);

}

void Concentrations::CreateCnE(mfem::ParGridFunction &Cn, double initial_value) {

    Cn = initial_value;

}

void Concentrations::Lithiation(mfem::ParGridFunction &Cn, double initial_value, ParFiniteElementSpace *fespace) {

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

void Concentrations::SBM_Matrix(mfem::ParGridFunction &psx, HypreParMatrix &Mmat, ParFiniteElementSpace *fespace) {

    // SetupBoundaryConditions();

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

void Concentrations::SetupBoundaryConditions(ParFiniteElementSpace *fespace) {
    
    Array<int> boundary_dofs;
    
    // Boundary attributes for Neumann BC on the west boundary
    // Array<int> nbc_w_bdr(pmesh->bdr_attributes.Max());
    nbc_w_bdr = 0; 
    nbc_w_bdr[0] = 1;  // Applying Neumann BC to the west boundary

    // Dirichlet BC on the east boundary for CnP
    Array<int> dbc_e_bdr(pmesh->bdr_attributes.Max());
    dbc_e_bdr = 0; 
    dbc_e_bdr[2] = 1;  // Applying Dirichlet BC to the east boundary

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

void Concentrations::ImposeNeumannBC(mfem::ParGridFunction &PGF, mfem::ParGridFunction &psx) {

    PGF = psx;
    PGF.Neg();

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

void Concentrations::ForceTerm(ParFiniteElementSpace *fespace, GridFunctionCoefficient cXx, mfem::ParLinearForm &Fxx, Array<int> boundary, ConstantCoefficient m, bool apply_boundary_conditions) {

    EnsureValidBoundaryAndFESpace();
    DebugBoundaryArray(boundary);

    std::unique_ptr<ParLinearForm> Bx2(new ParLinearForm(fespace));	

    Bx2->AddDomainIntegrator(new DomainLFIntegrator(cXx));
    // ConstantCoefficient temp_coef(1.0);

    if (apply_boundary_conditions) {

        Bx2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m), boundary);

    }

    Bx2->Assemble();

    Fxx = std::move(*Bx2);

    // Bx2->Print(std::cout);


}

void Concentrations::DebugBoundaryArray(const Array<int> &boundary)
{
    std::cout << "Boundary array contents: ";
    for (int i = 0; i < boundary.Size(); i++)
    {
        std::cout << boundary[i] << " ";
    }
    std::cout << std::endl;
}

void Concentrations::EnsureValidBoundaryAndFESpace()
{
    // Check if the mesh has boundary attributes defined
    if (!pmesh->bdr_attributes.Size())
    {
        std::cerr << "Error: The mesh does not have any boundary attributes defined." << std::endl;
        return;
    }

    std::cout << "Mesh boundary attributes size: " << pmesh->bdr_attributes.Size() << std::endl;
    std::cout << "Boundary attributes: ";
    for (int i = 0; i < pmesh->bdr_attributes.Size(); i++)
    {
        std::cout << pmesh->bdr_attributes[i] << " ";
    }
    std::cout << std::endl;

    // Now check that the finite element space (fespace) is valid and initialized with the mesh
    if (fespace->GetMesh() == nullptr)
    {
        std::cerr << "Error: The finite element space is not associated with a mesh." << std::endl;
        return;
    }

    std::cout << "Finite element space associated with the mesh is valid." << std::endl;
}



void Concentrations::TotalReaction(mfem::ParGridFunction &Rx, double xCrnt) {

    xCrnt = 0.0;
    Array<double> VtxVal(nC);
    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Rx.GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        xCrnt += EAvg(ei) * EVol(ei);
    }

    double geCrnt;
    MPI_Allreduce(&xCrnt, &geCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    infx = geCrnt / (mesh_handler.L_w);

    // cout << "infx: " << infx << std::endl;

}


// void Concentrations::TestFESpace(std::shared_ptr<ParFiniteElementSpace> fespace) {

//     ParGridFunction test_f(fespace.get());

//     // Check the value of test_f
//     std::cout << "Concentrations - First value of test_f: " << test_f(0) << std::endl;
// }
