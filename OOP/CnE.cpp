#include "CnE.hpp"
#include "Constants.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

CnE::CnE(MeshHandler &mesh_handler)
    : mesh_handler(mesh_handler),
      fespace(mesh_handler.GetFESpace()), 
      pse(fespace), // Initialize psi using fespace
      AvP(fespace),
      gtPse(mesh_handler.GetTotalPse()), 
      rho(Constants::rho), 
      Cr(Constants::Cr), 
      CnEGridFunction(fespace),
      L_w(mesh_handler.GetLw()) 
{
    // Initialize psi with the values from mesh_handler.GetPsi()
    pse = *mesh_handler.GetPse();
    AvP = *mesh_handler.GetAvP();
}

void CnE::Initialize() {
    double Ce0 = 0.001; // initial value
    CnEGridFunction = Ce0;

    // Degree of lithiation
    // double Xfr = 0.0;

	// // for imposing Neumann BC
    ParGridFunction PeR(fespace);
    PeR = pse;
    PeR.Neg();
    GridFunctionCoefficient matCoef_R(PeR);

    // SBM mass matrix
    //HypreParMatrix Mmatp;
    std::unique_ptr<ParBilinearForm> Me(new ParBilinearForm(fespace));
    GridFunctionCoefficient cPe(&pse);
    Me->AddDomainIntegrator(new MassIntegrator(cPe));
    Me->Assemble();
    Array<int> boundary_dofs;
    Me->FormSystemMatrix(boundary_dofs, Mmate);

    HypreParMatrix *TmatL, *TmatR;			// matrices for CN scheme


    HypreSmoother Me_prec;
    CGSolver Me_solver(MPI_COMM_WORLD);

    Me_solver.iterative_mode = false;
    Me_solver.SetRelTol(1e-7);
    Me_solver.SetAbsTol(0);
    Me_solver.SetMaxIter(200);
    Me_solver.SetPrintLevel(0);
    Me_prec.SetType(HypreSmoother::Jacobi);
    Me_solver.SetPreconditioner(Me_prec);

    std::cout << "Mmate rows after initialization: " << Mmate.NumRows() << ", cols: " << Mmate.NumCols() << std::endl;
   // Mmatp.Print("Simplified Mmatp");

    //HypreParMatrix *Tmatp;

    // SBM stiffness matrix
    ParGridFunction De(fespace);

    // HypreParMatrix Kmatp;

    // Force vector
    ParBilinearForm *Ke2;
    ParGridFunction Rxe(fespace);
    ParLinearForm *Be2;
    ParLinearForm Fet(fespace);
    //HypreParVector Fcb(fespace);

    // Create a Vector for CnP
    HypreParVector CeV0(fespace), CeVn(fespace), RHSe(fespace);

    // int nDof = CpV0.Size();

    // // Vector of psi
    // HypreParVector PsVc(fespace);
    // psi.GetTrueDofs(PsVc);

    // cout << "Degree of Lithiation: " << Xfr << std::endl;
}

void CnE::TimeStep(double dt) {
    
    // std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(fespace));
    // GridFunctionCoefficient cPs(&psi);
    // Mt->AddDomainIntegrator(new MassIntegrator(cPs));
    // Mt->Assemble();
    // Array<int> boundary_dofs;
    // Mt->FormSystemMatrix(boundary_dofs, Mmatp);
    
    
    HypreParVector Feb(fespace);

    //std::cout << "DEBUG: Before accessing Mmatp in TimeStep" << std::endl;
    //Mmatp.Print("TimeStep Mmatp");
    //std::cout << "DEBUG: After accessing Mmatp in TimeStep" << std::endl;

    ParGridFunction Rxe(fespace);
    Rxe = Rxn;
    Rxe *= (-1.0*Constants::t_minus);

    HypreParVector X1v(fespace);

    // int nV = fespace->GetNV(); // Get number of vertices

    // ParGridFunction Rxc(fespace);
    // Rxc = Rxn;
    // Rxc /= rho;
    GridFunctionCoefficient cAe(&Rxe);

    int nV = fespace->GetNV(); // Get number of vertices

    int nE = fespace->GetNE();
    int nC = pow(2, fespace->GetMesh()->Dimension());
    
    Vector EVol(nE);


    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = fespace->GetMesh()->GetElementVolume(ei);
    }

    Array<double> VtxVal(nC) ;				// values of corner vertices
	Vector EAvg(nE);						// element average	
    double val;

    // total reaction
	double eCrnt = 0.0;
    double geCrnt = 0.0;
    double infx = 0.0;

	ParGridFunction CeT(fespace);
	double CeC = 0.0;
	double CeAvg = 0.0;
	double gCeC = 0.0;	



	for (int ei = 0; ei < nE; ei++){	  
		Rxe.GetNodalValues(ei,VtxVal) ;
		val = 0.0;
		for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
		EAvg(ei) = val/nC;		 		  
		eCrnt += EAvg(ei)*EVol(ei);	
	} 	
	MPI_Allreduce(&eCrnt, &geCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	infx = geCrnt/L_w;

    // Neumann BC
	ConstantCoefficient nbcCoef(infx);
	ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    std::unique_ptr<ParLinearForm> Be2(new ParLinearForm(fespace));
    Be2->AddDomainIntegrator(new DomainLFIntegrator(cAe));
    Be2->Assemble();
    ParLinearForm Fet(fespace);
    Fet = std::move(*Be2);

    ParGridFunction De(fespace);
    // salt diffusivity in the electrolyte					
    for (int vi = 0; vi < nV; vi++){
        // appendix equation A-21
        De(vi) = pse(vi) * Constants::D0 * exp(-7.02-830 * CnEGridFunction(vi) + 50000 * CnEGridFunction(vi) * CnEGridFunction(vi));
    }
    GridFunctionCoefficient cDe(&De) ;

    std::unique_ptr<ParBilinearForm> Ke2(new ParBilinearForm(fespace));
    //HypreParMatrix Kmate;

    Ke2->AddDomainIntegrator(new DiffusionIntegrator(cDe));
    Ke2->Assemble();
    Ke2->FormLinearSystem(boundary_dofs, CnEGridFunction, Fet, Kmate, X1v, Feb);
    Feb *= dt;

    //HypreParMatrix *Tmatp;

    // std::cout << "Boundary DOFs: " << boundary_dofs.Size() << endl;
    std::cout << "Mmate rows: " << Mmate.NumRows() << ", cols: " << Mmate.NumCols() << std::endl;
    std::cout << "Kmate rows: " << Kmate.NumRows() << ", cols: " << Kmate.NumCols() << std::endl;

    // Mmatp.Print("Mmatp");
    // Kmatp.Print("Kmatp");

// // Crank-Nicolson matrices	
    TmatR = Add(1.0, Mmate, -0.5*dt, Kmate);		
    TmatL = Add(1.0, Mmate,  0.5*dt, Kmate);
    
    // vector of CnE				
    CnE.GetTrueDofs(CeV0);		
            
    TmatR->Mult(CeV0,RHSe);
    RHSe += Feb;
    
    // solver
    Me_solver.SetOperator(*TmatL);    	
    
    // time stepping
    Me_solver.Mult(RHSe,CeVn) ;
    
    // recover
    CnE.Distribute(CeVn);    	

    // check conservation of salt
    if (t%500 == 0 && t > 0){
        CeC = 0.0;
        CeT = CnE;
        CeT *= pse;
        for (int ei = 0; ei < nE; ei++){
            CeT.GetNodalValues(ei,VtxVal) ;
            val = 0.0;
            for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
            EAvg(ei) = val/nC;	
            CeC += EAvg(ei)*EVol(ei) ;
        }
        MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
        // average CnE throughout electrolyte
        CeAvg = gCeC/gtPse;				
        
        // adjust CnE
        CnE -= (CeAvg-Ce0);
        MPI_Barrier(MPI_COMM_WORLD);
    }	
    TmatR->~HypreParMatrix();
    TmatL->~HypreParMatrix();		
    Be2->~ParLinearForm();

    std::cout << "Tmatp rows: " << Tmatp->NumRows() << ", cols: " << Tmatp->NumCols() << std::endl;
}

void CnE::Save() {
    CnEGridFunction.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/OOPCnE");
}


