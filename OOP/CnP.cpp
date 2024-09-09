#include "CnP.hpp"
#include "MeshHandler.hpp"
#include "Constants.hpp"
#include <fstream>
#include <iostream>
#include "PotP.hpp"

using namespace mfem;
using namespace std;

CnP::CnP(MeshHandler &mesh_handler, PotP &potp)
    : mesh_handler(mesh_handler),
      pmesh(mesh_handler.GetParMesh()),
      fespace(mesh_handler.GetFESpace()), 
      psi(fespace), 
      AvP(fespace),
      AvB(fespace),
      gtPsi(mesh_handler.GetTotalPsi()), 
      rho(Constants::rho), 
      Cr(Constants::Cr), 
      CnPGridFunction(fespace),
      TmpF(fespace),
      potp(potp),
      kap(fespace),
      B1t(fespace),
      B1v(fespace),
      phP(potp.GetphP()) // defined in PotP.hpp

{
    // Initialize psi with the values from mesh_handler.GetPsi()
    psi = *mesh_handler.GetPsi();
    AvP = *mesh_handler.GetAvP();
    AvB = *mesh_handler.GetAvB();

    BvP = potp.GetBvP();
}

void CnP::Initialize() {

    Cp0 = 0.3; // initial value
    CnPGridFunction = Cp0;

    // Degree of lithiation
    Xfr = 0.0;

    TmpF = CnPGridFunction;
    TmpF *= psi;

    lSum = 0.0;
    nE = fespace->GetNE(); // Get number of elements
    nC = pow(2, fespace->GetMesh()->Dimension()); // Number of corner vertices

    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = fespace->GetMesh()->GetElementVolume(ei);
    }
    Vector EAvg(nE);
    Array<double> VtxVal(nC);

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

    MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Xfr = gSum / gtPsi;

    // SBM mass matrix
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

    // SBM stiffness matrix
    ParGridFunction Dp(fespace);

    // Force vector
    ParGridFunction Rxc(fespace);
    ParLinearForm *Bc2;
    ParLinearForm Fct(fespace);

    // Create a Vector for CnP
    HypreParVector CpV0(fespace), CpVn(fespace), RHCp(fespace);

    nDof = CpV0.Size();

    // Vector of psi
    HypreParVector PsVc(fespace);
    psi.GetTrueDofs(PsVc);

    cout << "Degree of Lithiation: " << Xfr << std::endl;

    Rxn = make_unique<ParGridFunction>(fespace); // this was the difference!
    *Rxn = AvP; 
    *Rxn *= 1.0e-10; 

}

void CnP::TimeStep(double dt) {
    
    // make sub functions to clean this portion up
    std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(fespace));
    GridFunctionCoefficient cPs(&psi);
    Mt->AddDomainIntegrator(new MassIntegrator(cPs));
    Mt->Assemble();
    Array<int> boundary_dofs;
    Mt->FormSystemMatrix(boundary_dofs, Mmatp);
      
    HypreParVector Fcb(fespace);

    ParGridFunction Rxn(fespace);

    HypreParVector X1v(fespace);

    nV = fespace->GetNV(); // Get number of vertices

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

    // try just doing Kc2 update instead of lines 150-152
    
    std::unique_ptr<ParBilinearForm> Kc2(new ParBilinearForm(fespace));

    Kc2->AddDomainIntegrator(new DiffusionIntegrator(cDp));
    Kc2->Assemble();
    Kc2->FormLinearSystem(boundary_dofs, CnPGridFunction, Fct, Kmatp, X1v, Fcb);
    Fcb *= dt;

    std::cout << "Boundary DOFs: " << boundary_dofs.Size() << endl;
    std::cout << "Mmatp rows: " << Mmatp.NumRows() << ", cols: " << Mmatp.NumCols() << std::endl;
    std::cout << "Kmatp rows: " << Kmatp.NumRows() << ", cols: " << Kmatp.NumCols() << std::endl;

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

    lSum = 0.0;
    nE = fespace->GetNE();
    nC = pow(2, fespace->GetMesh()->Dimension());
    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = fespace->GetMesh()->GetElementVolume(ei);
    }
    Vector EAvg(nE);
    Array<double> VtxVal(nC);

    for (int ei = 0; ei < nE; ei++) {
        TmpF.GetNodalValues(ei, VtxVal);
        val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        lSum += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Xfr = gSum / gtPsi;

    delete Tmatp;

    // STARTING FROM REACTION PART

    // assign known values to the DBC nodes	
    ConstantCoefficient dbc_e_Coef(BvP);			
    
    // particle conductivity
    // appendix equation A-20
    for (int vi = 0; vi < nV; vi++){
        kap(vi) = psi(vi)*(0.01929 + 0.7045*tanh(2.399*CnPGridFunction(vi)) - \
            0.7238*tanh(2.412*CnPGridFunction(vi)) - 4.2106e-6);
    }	
    GridFunctionCoefficient cKp(&kap) ;
		
    // stiffness matrix for phP
    std::unique_ptr<ParBilinearForm> Kp2(new ParBilinearForm(fespace));

    Kp2->AddDomainIntegrator(new DiffusionIntegrator(cKp));
    Kp2->Assemble();

    // Dirichlet BC on the east boundary. phP
	Array<int> dbc_e_bdr(pmesh->bdr_attributes.Max());
	dbc_e_bdr = 0; dbc_e_bdr[2] = 1;
	// use dbc_e_bdr array to extract all node labels of Dirichlet BC
	Array<int> ess_tdof_list_e(0);			
	fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);
    
    // project values to DBC nodes
    phP.ProjectBdrCoefficient(dbc_e_Coef, dbc_e_bdr); 	
    Kp2->FormLinearSystem(ess_tdof_list_e, phP, B1t, KmP, X1v, B1v);			

    // Solve the system using PCG with hypre's BoomerAMG preconditioner.
    
    CGSolver cgPP(MPI_COMM_WORLD);
	cgPP.SetRelTol(1e-7);
	cgPP.SetMaxIter(200);
    
    HypreBoomerAMG Mpp(KmP);
    Mpp.SetPrintLevel(0);
    cgPP.SetPreconditioner(Mpp);
    cgPP.SetOperator(KmP);

    // rate constants and exchange current density at interface
    for (int vi = 0; vi < nV; vi++){
        if ( AvB(vi)*Constants::dh > 0.0 ){
            val = -0.2*(CnPGridFunction(vi)-0.37)-1.559-0.9376*tanh(8.961*CnPGridFunction(vi)-3.195);
            i0C(vi) = pow(10.0,val)*1.0e-3;
            
            OCV(vi) = 1.095*CnP(vi)*CnPGridFunction(vi) - 8.324e-7*exp(14.31*CnPGridFunction(vi)) + \
                4.692*exp(-0.5389*CnPGridFunction(vi));
                
            Kfw(vi) = i0C(vi)/(Constants::Frd*0.001  )*exp( Constants::alp*Constants::Cst1*OCV(vi)) ;	
            Kbw(vi) = i0C(vi)/(Constants::Frd*CnPGridFunction(vi))*exp(-Constants::alp*Constants::Cst1*OCV(vi)) ;
        }
    }



}

void CnP::Save() {
    CnPGridFunction.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/OOPCnP");
}



