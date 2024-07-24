#include "CnP.hpp"
#include "Constants.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

CnP::CnP(MeshHandler &mesh_handler)
    : mesh_handler(mesh_handler),
      fespace(mesh_handler.GetFESpace()), 
      psi(fespace), // Initialize psi using fespace
      AvP(fespace),
      gtPsi(mesh_handler.GetTotalPsi()), 
      rho(Constants::rho), 
      Cr(Constants::Cr), 
      CnPGridFunction(fespace),
      Rxn(std::make_unique<mfem::ParGridFunction>(fespace)) // Initialize Rxn

{
    // Initialize psi with the values from mesh_handler.GetPsi()
    psi = *mesh_handler.GetPsi();
    AvP = *mesh_handler.GetAvP();
}

void CnP::Initialize() {
    double Cp0 = 0.3; // initial value
    CnPGridFunction = Cp0;

    // Degree of lithiation
    double Xfr = 0.0;

    ParGridFunction TmpF(fespace);
    TmpF = CnPGridFunction;
    TmpF *= psi;

    double lSum = 0.0;
    int nE = fespace->GetNE(); // Get number of elements
    int nC = pow(2, fespace->GetMesh()->Dimension()); // Number of corner vertices
    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = fespace->GetMesh()->GetElementVolume(ei);
    }
    Vector EAvg(nE);
    Array<double> VtxVal(nC);
    double val;

    for (int ei = 0; ei < nE; ei++) {
        TmpF.GetNodalValues(ei, VtxVal);
        val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        lSum += EAvg(ei) * EVol(ei);
    }

    cout << "Sum: " << lSum << std::endl;

    double gSum;
    MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Xfr = gSum / gtPsi;

    // SBM mass matrix
    //HypreParMatrix Mmatp;
    std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(fespace));
    GridFunctionCoefficient cPs(&psi);
    Mt->AddDomainIntegrator(new MassIntegrator(cPs));
    Mt->Assemble();
    Array<int> boundary_dofs;
    Mt->FormSystemMatrix(boundary_dofs, Mmatp);

    HypreSmoother Mp_prec;
    CGSolver Mp_solver(MPI_COMM_WORLD);

    Mp_solver.iterative_mode = false;
    Mp_solver.SetRelTol(1e-7);
    Mp_solver.SetAbsTol(0);
    Mp_solver.SetMaxIter(200);
    Mp_solver.SetPrintLevel(0);
    Mp_prec.SetType(HypreSmoother::Jacobi);
    Mp_solver.SetPreconditioner(Mp_prec);
    Mp_solver.SetOperator(Mmatp);

    std::cout << "Mmatp rows after initialization: " << Mmatp.NumRows() << ", cols: " << Mmatp.NumCols() << std::endl;
   // Mmatp.Print("Simplified Mmatp");

    //HypreParMatrix *Tmatp;

    // SBM stiffness matrix
    ParGridFunction Dp(fespace);

    // HypreParMatrix Kmatp;

    // Force vector
    ParGridFunction Rxc(fespace);
    ParLinearForm *Bc2;
    ParLinearForm Fct(fespace);
    //HypreParVector Fcb(fespace);

    // Create a Vector for CnP
    HypreParVector CpV0(fespace), CpVn(fespace), RHCp(fespace);

    int nDof = CpV0.Size();

    // Vector of psi
    HypreParVector PsVc(fespace);
    psi.GetTrueDofs(PsVc);

    cout << "Degree of Lithiation: " << Xfr << std::endl;
}

void CnP::TimeStep(double dt) {
    
    std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(fespace));
    GridFunctionCoefficient cPs(&psi);
    Mt->AddDomainIntegrator(new MassIntegrator(cPs));
    Mt->Assemble();
    Array<int> boundary_dofs;
    Mt->FormSystemMatrix(boundary_dofs, Mmatp);
    
    
    HypreParVector Fcb(fespace);

    //std::cout << "DEBUG: Before accessing Mmatp in TimeStep" << std::endl;
    //Mmatp.Print("TimeStep Mmatp");
    //std::cout << "DEBUG: After accessing Mmatp in TimeStep" << std::endl;

    ParGridFunction Rxn(fespace);
    Rxn = 0.0;
    Rxn = AvP;
    Rxn *= 1.0e-10;

    HypreParVector X1v(fespace);

    int nV = fespace->GetNV(); // Get number of vertices

    ParGridFunction Rxc(fespace);
    Rxc = *Rxn;
    Rxc /= rho;
    GridFunctionCoefficient cAp(&Rxc);

    std::unique_ptr<ParLinearForm> Bc2(new ParLinearForm(fespace));
    Bc2->AddDomainIntegrator(new DomainLFIntegrator(cAp));
    Bc2->Assemble();
    ParLinearForm Fct(fespace);
    Fct = std::move(*Bc2);

    ParGridFunction Dp(fespace);
    for (int vi = 0; vi < nV; vi++) {
        Dp(vi) = psi(vi) * (0.0277 - 0.084 * CnPGridFunction(vi) + 0.1003 * CnPGridFunction(vi) * CnPGridFunction(vi)) * 1.0e-8;
        if (Dp(vi) > 4.6e-10) {
            Dp(vi) = 4.6e-10;
        }
    }
    GridFunctionCoefficient cDp(&Dp);

    std::unique_ptr<ParBilinearForm> Kc2(new ParBilinearForm(fespace));
    //HypreParMatrix Kmatp;

    Kc2->AddDomainIntegrator(new DiffusionIntegrator(cDp));
    Kc2->Assemble();
    Kc2->FormLinearSystem(boundary_dofs, CnPGridFunction, Fct, Kmatp, X1v, Fcb);
    Fcb *= dt;

    //HypreParMatrix *Tmatp;

    std::cout << "Boundary DOFs: " << boundary_dofs.Size() << endl;
    std::cout << "Mmatp rows: " << Mmatp.NumRows() << ", cols: " << Mmatp.NumCols() << std::endl;
    std::cout << "Kmatp rows: " << Kmatp.NumRows() << ", cols: " << Kmatp.NumCols() << std::endl;

    // Mmatp.Print("Mmatp");
    // Kmatp.Print("Kmatp");

    Tmatp = Add(1.0, Mmatp, -dt, Kmatp);



    std::cout << "Tmatp rows: " << Tmatp->NumRows() << ", cols: " << Tmatp->NumCols() << std::endl;

    HypreParVector CpV0(fespace), CpVn(fespace), RHCp(fespace);
    CnPGridFunction.GetTrueDofs(CpV0);

    Tmatp->Mult(CpV0, RHCp);
    RHCp += Fcb;


    HypreSmoother Mp_prec;
    CGSolver Mp_solver(MPI_COMM_WORLD);

    Mp_solver.iterative_mode = false;
    Mp_solver.SetRelTol(1e-7);
    Mp_solver.SetAbsTol(0);
    Mp_solver.SetMaxIter(200);
    Mp_solver.SetPrintLevel(0);
    Mp_prec.SetType(HypreSmoother::Jacobi);
    Mp_solver.SetPreconditioner(Mp_prec);
    Mp_solver.SetOperator(Mmatp);

    Mp_solver.Mult(RHCp, CpVn);

    HypreParVector PsVc(fespace);
    psi.GetTrueDofs(PsVc);

    int nDof = CpV0.Size();
    for (int p = 0; p < nDof; p++) {
        if (PsVc(p) < 1.0e-5) {
            CpVn(p) = 0.3; // Cp0
        }
    }

    CnPGridFunction.Distribute(CpVn);

    ParGridFunction TmpF(fespace);
    TmpF = CnPGridFunction;
    TmpF *= psi;

    double lSum = 0.0;
    int nE = fespace->GetNE();
    int nC = pow(2, fespace->GetMesh()->Dimension());
    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = fespace->GetMesh()->GetElementVolume(ei);
    }
    Vector EAvg(nE);
    Array<double> VtxVal(nC);
    double val;

    for (int ei = 0; ei < nE; ei++) {
        TmpF.GetNodalValues(ei, VtxVal);
        val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        lSum += EAvg(ei) * EVol(ei);
    }

    double gSum;
    MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double Xfr = gSum / gtPsi;
    
    delete Tmatp;

    // cout << "Degree of Lithiation: " << Xfr << std::endl;
}

void CnP::Save() {
    CnPGridFunction.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/OOPCnP");
}



