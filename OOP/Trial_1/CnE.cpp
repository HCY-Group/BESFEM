#include "CnE.hpp"
#include "CnP.hpp"
#include "PotE.hpp"
#include "Constants.hpp"
#include "Reaction.hpp"
#include <iostream>
#include "MeshHandler.hpp" 
// #include "mfem.hpp"

using namespace mfem;
using namespace std;

CnE::CnE(MeshHandler &mesh_handler, Reaction &reaction, CnP &cnp)
    : mesh_handler(mesh_handler),
      pmesh(mesh_handler.GetParMesh()),
      fespace(mesh_handler.GetFESpace()), 
      reaction(reaction),
      pse(fespace), // Initialize psi using fespace
      AvP(fespace),
      gtPse(mesh_handler.GetTotalPse()), 
      rho(Constants::rho), 
      Cr(Constants::Cr), 
      CnEGridFunction(fespace),
      L_w(mesh_handler.GetLw()),
      cnp(cnp),
      Dmp(fespace),
      kpl(fespace),
      CeVn(fespace),
      B1t(fespace),
      X1v(fespace),
      B1v(fespace),
      LpCe(fespace)

{

    // Initialize psi with the values from mesh_handler.GetPsi()
    pse = *mesh_handler.GetPse();
    AvP = *mesh_handler.GetAvP();

    // tc1 = pote.Gettc1();
    // tc2 = pote.Gettc2();

    Rxn = *cnp.GetRxn();

    BvE = -1.0;

}


void CnE::Initialize() {
    // std::cout << "Entering CnE::Initialize()" << std::endl;

    double Ce0 = 0.001; // initial value
    CnEGridFunction = Ce0;

	// // for imposing Neumann BC
    ParGridFunction PeR(fespace);
    PeR = pse;
    PeR.Neg();
    GridFunctionCoefficient matCoef_R(&PeR);


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

    // SBM stiffness matrix
    ParGridFunction De(fespace);

    // Force vector
    ParBilinearForm *Ke2;
    ParGridFunction Rxe(fespace);
    ParLinearForm *Be2;
    ParLinearForm Fet(fespace);

    int nE = fespace->GetNE(); // Get number of elements
    int nC = pow(2, fespace->GetMesh()->Dimension()); // Number of corner vertices

    Array<double> VtxVal(nC);
    Vector EVol(nE);
    Vector EAvg(nE);
    double eCrnt = 0.0;

    double geCrnt = 0.0;
	double infx = 0.0;

    double CeC = 0.0;
    double gCeC = 0.0;
    double CeAvg = 0.0;

    ParGridFunction CeT(fespace);

    int s = 0;

}

void CnE::TimeStep(double dt) {

    std::unique_ptr<ParBilinearForm> Me(new ParBilinearForm(fespace));
    GridFunctionCoefficient cPe(&pse);
    Me->AddDomainIntegrator(new MassIntegrator(cPe));
    Me->Assemble();
    Array<int> boundary_dofs;
    Me->FormSystemMatrix(boundary_dofs, Mmate);
    
    HypreParVector Feb(fespace);
    HypreParVector CeV0(fespace), RHSe(fespace);

    // cout << "Rxn bringing into CnE" << endl;
    // Rxn.Print(std::cout);

    // ParGridFunction Rxe(fespace);
    Rxe = make_unique<ParGridFunction>(fespace);
    *Rxe = Rxn;
    *Rxe *= (-1.0*Constants::t_minus);

    // GridFunctionCoefficient cAe(&Rxe);
    GridFunctionCoefficient cAe(Rxe.get()); // had to pass pointer from unique_ptr to GFC

    nC = pow(2, fespace->GetMesh()->Dimension()); // Number of corner vertices

    nE = fespace->GetNE(); // Get number of elements

    nV = fespace->GetNV();

    Array<double> VtxVal(nC);
    Vector EVol(nE);
	for (int ei = 0; ei < nE; ei++){
		// EVol(ei) = pmesh.GetElementVolume(ei);	
        EVol(ei) = fespace->GetMesh()->GetElementVolume(ei);

	} 	
    Vector EAvg(nE);

	eCrnt = 0.0;
    for (int ei = 0; ei < nE; ei++){

		// Rxe.GetNodalValues(ei,VtxVal) ;
        Rxe->GetNodalValues(ei, VtxVal);  // -> due to unique_ptr
		val = 0.0;
		for (int vt = 0; vt < nC; vt++){
            val += VtxVal[vt];
        }
		EAvg(ei) = val/nC;		 		  
		eCrnt += EAvg(ei)*EVol(ei);	
	} 	

	MPI_Allreduce(&eCrnt, &geCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	infx = geCrnt/L_w;

    cout << "eCrnt: " << eCrnt << endl;

    // Neumann BC
	ConstantCoefficient nbcCoef(infx);
	ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

 	// force term
    std::unique_ptr<ParLinearForm> Be2(new ParLinearForm(fespace));
    Be2->AddDomainIntegrator(new DomainLFIntegrator(cAe));
    // Be2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m_nbcCoef), nbc_w_bdr);
    Be2->Assemble();
    ParLinearForm Fet(fespace);
    Fet = std::move(*Be2);
    // Rxe.Print(std::cout);

    ParGridFunction De(fespace);
    // salt diffusivity in the electrolyte					
    for (int vi = 0; vi < nV; vi++){
        // appendix equation A-21
        De(vi) = pse(vi) * Constants::D0 * exp(-7.02-830 * CnEGridFunction(vi) + 50000 * CnEGridFunction(vi) * CnEGridFunction(vi));
    }
    GridFunctionCoefficient cDe(&De) ;

    // K Matrix

    std::unique_ptr<ParBilinearForm> Ke2(new ParBilinearForm(fespace));

    HypreParVector X1v(fespace);

    Ke2->AddDomainIntegrator(new DiffusionIntegrator(cDe));
    Ke2->Assemble();
    Ke2->FormLinearSystem(boundary_dofs, CnEGridFunction, Fet, Kmate, X1v, Feb);
    Feb *= dt;

    // Crank-Nicolson matrices	
    TmatR = Add(1.0, Mmate, -0.5*dt, Kmate);		
    TmatL = Add(1.0, Mmate,  0.5*dt, Kmate);
    
    // vector of CnE				
    CnEGridFunction.GetTrueDofs(CeV0);		
            
    TmatR->Mult(CeV0,RHSe);

    // solver

    HypreSmoother Me_prec;
    CGSolver Me_solver(MPI_COMM_WORLD);

    Me_solver.iterative_mode = false;
    Me_solver.SetRelTol(1e-7);
    Me_solver.SetAbsTol(0);
    Me_solver.SetMaxIter(200);
    Me_solver.SetPrintLevel(0);
    Me_prec.SetType(HypreSmoother::Jacobi);
    Me_solver.SetPreconditioner(Me_prec);
    
    Me_solver.SetOperator(*TmatL);    	
    Me_solver.Mult(RHSe, CeVn); 

    // recover
    CnEGridFunction.Distribute(CeVn);   
    ParGridFunction CeT(fespace);

    Ce0 = 0.001;
 	    
    // check conservation of salt
    if (s%1 == 0 && s > 0){ // used variable s for timestep 
        CeC = 0.0;
        CeT = CnEGridFunction;
        CeT *= pse;
        for (int ei = 0; ei < nE; ei++){
            CeT.GetNodalValues(ei,VtxVal) ;
            val = 0.0;
            for (int vt = 0; vt < nC; vt++){
                val += VtxVal[vt];
            }
            EAvg(ei) = val/nC;	
            CeC += EAvg(ei)*EVol(ei) ;
        }
        
        cout << "CeC: " << CeC << endl;

        MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
        // average CnE throughout electrolyte
        CeAvg = gCeC/gtPse;	
 
        // adjust CnE
        CnEGridFunction -= (CeAvg-Ce0);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    cout << "CeAvg: " << CeAvg << endl;

    s = s + 1;



    // // STARTING FROM REACTION PART

    tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
	tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    cout << "tc1 :" << tc1 << endl;

    // electrolyte conductivity and RHS	
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02-830*CnEGridFunction(vi)+50000*CnEGridFunction(vi)*CnEGridFunction(vi));
        Dmp(vi) = pse(vi)*tc1*Constants::D0*dffe;
        kpl(vi) = pse(vi)*tc2*Constants::D0*dffe*CnEGridFunction(vi);
    }
    GridFunctionCoefficient cDm(&Dmp);

    // Dmp.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/DmpOOP_CnE");

}

void CnE::Save() {
    CnEGridFunction.Save("/mnt/home/brandlan/PhD/MFEM_Parallel/mfem-4.5/GitLab/besfem/OOP/OOPCnE");
}


