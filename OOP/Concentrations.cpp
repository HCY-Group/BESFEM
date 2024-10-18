#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include <fstream>
#include <iostream>

#include <HYPRE.h>
#include <HYPRE_utilities.h>
#include "mfem.hpp"



using namespace mfem;
using namespace std;

Concentrations::Concentrations(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : mesh_handler(mh), fespace(fe), pmesh(pm), psi(*mh.GetPsi()), pse(*mh.GetPse()), EVol(mh.GetElementVolume()),
    reaction(fe, mh, *this), Mmat(nullptr), Mmatp(nullptr),
    Mp_solver(std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD))



    // , psi(*mesh_handler.GetPsi()), pse(*mesh_handler.GetPse()), , EVol(mesh_handler.GetElementVolume()),  M(nullptr), , solver(nullptr)
    //   EVol(mesh_handler.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi()), reaction(mesh_handler, *this), PeR(fe), matCoef_R(&PeR), Mp_solver(nullptr)

      
{

    // nbc_w_bdr.SetSize(mesh_handler.pmesh->bdr_attributes.Max());
    // nbc_w_bdr = 0;  // Initialize the array with zeros

    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();
    
    // Initialize ParGridFunction for CnP and CnE
    CnP = new ParGridFunction(fespace);
    CnE = new ParGridFunction(fespace);

    Rxc = new ParGridFunction(fespace);
    Rxe = new ParGridFunction(fespace);
    PeR = new ParGridFunction(fespace);
    // Rxe = make_unique<ParGridFunction>(fespace.get());

    // X1v = HypreParVector(fespace.get());
    Fcb = HypreParVector(fespace);
    PsVc = HypreParVector(fespace);
    // Feb = HypreParVector(fespace.get());

    // reaction.Initialize(); 

    // Mp_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    // Me_solver = mfem::CGSolver(MPI_COMM_WORLD);

}


void Concentrations::InitializeCnP() {

    Lithiation(*CnP, 0.3);

    Mmatp = std::make_shared<mfem::HypreParMatrix>();
    Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);

    Solver(psi, Mmatp, *Mp_solver, Mp_prec);
    
    // cout << "CnP in Initialize Step:" << endl;
    // CnP->Print(std::cout);

	psi.GetTrueDofs(PsVc);
    double *PsVc_data = PsVc.GetData();

}


void Concentrations::InitializeCnE() {

    CreateCnE(*CnE, 0.001);

    Mmate = std::make_shared<mfem::HypreParMatrix>();
    Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);

    Solver(pse, Mmate, *Me_solver, Me_prec);
    
    // ParGridFunction PeR(fespace);
	// PeR = pse;
	// PeR.Neg();
    // GridFunctionCoefficient matCoef_R(&PeR);
    ImposeNeumannBC(pse, *PeR);

    cout << "CnE in Initialize Step:" << endl;
    CnE ->Print(std::cout);

}


void Concentrations::TimeStepCnP() {
    
    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    mfem::GridFunctionCoefficient cAp(Rxc);
    SetupRx(Rxn, *Rxc, Constants::rho, cAp);

    Array<int> dummy_boundary;
    ConstantCoefficient coef1(0.0);
    ConstantCoefficient coef2(0.0);
    ProductCoefficient dummy_coef(coef1, coef2);

    ForceTerm(*Rxc, Fct, dummy_boundary, dummy_coef, false); // false since not applying BCs

    std::shared_ptr<GridFunctionCoefficient> cDp = Diffusivity(psi, *CnP, true); // true since using first equation

    Kmatp = std::make_shared<mfem::HypreParMatrix>();
    K_Matrix(boundary_dofs, *CnP, Fct, Kmatp, X1v, Fcb, cDp.get());

    // Create T Matrix
    Tmatp = Add(1.0, *Mmatp, -(Constants::dt), *Kmatp);

    HypreParVector CpV0(fespace);
    int nDof = CpV0.Size();
    CnP->GetTrueDofs(CpV0);

    HypreParVector RHCp(fespace);
    Tmatp->Mult(CpV0, RHCp);
    RHCp += Fcb;

    HypreParVector CpVn(fespace);

    Mp_solver->Mult(RHCp, CpVn);

    // Update only the solid region
    for (int p = 0; p < nDof; p++){
        if (PsVc(p) < 1.0e-5){
            CpVn(p) = 0.3;} // Cp0 initial value
    }

    CnP->Distribute(CpVn);

    // Degree of Lithiation
    LithiationCalculation(*CnP);

    // std::cout << "Updated CnP values:" << std::endl;
    // CnP->Print(std::cout);

}

void Concentrations::TimeStepCnE() {

    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    GridFunctionCoefficient cAe(Rxe);
    SetupRx(Rxn, *Rxe, Constants::t_minus, cAe);

    TotalReaction(*Rxe, eCrnt);
    ConstantCoefficient nbcCoef(infx); // Neumann BC works with -infx

    GridFunctionCoefficient matCoef_R(PeR);

    ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    Array<int> nbc_w_bdr(pmesh->bdr_attributes.Max());
	nbc_w_bdr = 0; 
	nbc_w_bdr[0] = 1;

    ForceTerm(*Rxe, Fet, nbc_w_bdr, m_nbcCoef, true); // true since applying boundary conditions

    std::shared_ptr<GridFunctionCoefficient> cDe = Diffusivity(pse, *CnE, false); // false using other equation
    
    Kmate = std::make_shared<mfem::HypreParMatrix>();
    K_Matrix(boundary_dofs, *CnE, Fet, Kmate, X1v, Feb, cDe.get());

    // Crank-Nicolson matrices
    TmatR = Add(1.0,*Mmate, -0.5*Constants::dt, *Kmate);		
    TmatL = Add(1.0, *Mmate,  0.5*Constants::dt, *Kmate);	
    
    HypreParVector CeV0(fespace);
    CnE->GetTrueDofs(CeV0);		

    HypreParVector RHCe(fespace);
    TmatR->Mult(CeV0, RHCe);
    RHCe += Feb;
    
    HypreParVector CeVn(fespace);

    Me_solver->SetOperator(*TmatL);
    Me_solver->Mult(RHCe, CeVn) ;

	CnE->Distribute(CeVn);    	


    std::cout << "Updated CnE values:" << std::endl;
    CnE->Print(std::cout);

}

void Concentrations::CreateCnE(mfem::ParGridFunction &Cn, double initial_value) {

    for (int i = 0; i < Cn.Size(); ++i) {
        Cn(i) = initial_value;  // Set all values of Cn to initial_value
    }    

}

void Concentrations::LithiationCalculation(mfem::ParGridFunction &Cn) {

    // Cn = initial_value;

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

void Concentrations::Lithiation(mfem::ParGridFunction &Cn, double initial_value) {

    for (int i = 0; i < Cn.Size(); ++i) {
        Cn(i) = initial_value;  // Set all values of Cn to initial_value
    }    
    LithiationCalculation(Cn);
}


// want to use the Mmatp here as was used in SBM function
void Concentrations::Solver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother) {
    
    M = new ParBilinearForm(fespace);

    Ps_gf = new ParGridFunction(fespace);
    *Ps_gf = psx;

    cP = new GridFunctionCoefficient(Ps_gf);

    M->AddDomainIntegrator(new MassIntegrator(*cP));
    M->Assemble();
    M->Finalize();

    HypreParMatrix HPM;
    M->FormSystemMatrix(boundary_dofs, HPM);
    Mmat = std::make_shared<mfem::HypreParMatrix>(HPM);

    // mfem::HypreSmoother M_prec;
    smoother.SetType(mfem::HypreSmoother::Jacobi);

    m_solver.iterative_mode = false;
    m_solver.SetRelTol(1e-7);
    m_solver.SetAbsTol(0);
    m_solver.SetMaxIter(102);
    m_solver.SetPrintLevel(0);
    m_solver.SetPreconditioner(smoother);
    m_solver.SetOperator(*Mmat);

}

void Concentrations::ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF) {

    PGF = psx;
    PGF.Neg();

}

// std::shared_ptr<mfem::GridFunctionCoefficient> Concentrations::ImposeNeumannBC(mfem::ParGridFunction &psx){

//     // mfem::ParGridFunction* coef = new mfem::ParGridFunction(fespace);

//     std::shared_ptr<mfem::ParGridFunction> PGF = std::make_shared<mfem::ParGridFunction>(fespace);
//     *PGF = psx;
//     PGF->Neg();

//     coef = std::make_shared<mfem::GridFunctionCoefficient>(PGF.get());


//     return coef;

// }

void Concentrations::SetupRx(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, 
double value, GridFunctionCoefficient cAx) {

    Rx2 = Rx1;

    if (&Rx2 == Rxc) {
        Rx2 /= value; // this is needed for CnP & Rxc
    }

    if (&Rx2 == Rxe) {
        Rx2 *= (-1.0 * value); // this is needed for CnE & Rxe
    }
    
}

void Concentrations::ForceTerm(mfem::ParGridFunction &gfc, mfem::ParLinearForm &Fxx, Array<int> boundary, ProductCoefficient m, bool apply_boundary_conditions) {

    std::unique_ptr<ParLinearForm> Bx2(new ParLinearForm(fespace));	
    // Bx2 = new ParLinearForm(fespace);

    Rxx = new ParGridFunction(fespace);
    *Rxx = gfc;

    cXx = new GridFunctionCoefficient(Rxx);
    Bx2->AddDomainIntegrator(new DomainLFIntegrator(*cXx));

    if (apply_boundary_conditions) {
        Bx2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m), boundary);
    }

    Bx2->Assemble();
    Fxx = std::move(*Bx2);

    // cout << "Bx2 in Concentrations" << endl; // fix this
    // Bx2->Print(std::cout);

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

}


std::shared_ptr<mfem::GridFunctionCoefficient> Concentrations::Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte ){

    mfem::ParGridFunction* Dx = new mfem::ParGridFunction(fespace);

    for (int vi = 0; vi < nV; vi++) {
        if (particle_electrolyte) {
            (*Dx)(vi) = psx(vi) * (0.0277 - 0.084 * Cn(vi) + 0.1003 * Cn(vi) * Cn(vi)) * 1.0e-8;
            if ((*Dx)(vi) > 4.6e-10) {
                (*Dx)(vi) = 4.6e-10;
            }
        } else {
            (*Dx)(vi) = psx(vi) * Constants::D0 * exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
        }
    }

    return std::make_shared<mfem::GridFunctionCoefficient>(Dx);

}

void Concentrations::K_Matrix(Array<int> boundary, mfem::ParGridFunction &Cn, ParLinearForm &Fxx, std::shared_ptr<HypreParMatrix> &Kmatx, HypreParVector &X1v, HypreParVector &Fxb, GridFunctionCoefficient *cDx) {

    // SetupBoundaryConditions();
    std::unique_ptr<ParBilinearForm> Kx2(new ParBilinearForm(fespace));

    HypreParMatrix Khpm;
    
    Kx2->AddDomainIntegrator(new DiffusionIntegrator(*cDx));
    Kx2->Assemble();
    Kx2->FormLinearSystem(boundary, Cn, Fxx, Khpm, X1v, Fxb);

    Kmatx = std::make_shared<mfem::HypreParMatrix>(Khpm);

    Fxb *= Constants::dt;

    // // Get the local data of the HypreParVector
    // double *Fxb_data = Fxb.GetData();

    // // Print each value of the vector
    // int size = Fxb.Size();
    // std::cout << "Fxb values:" << std::endl;
    // for (int i = 0; i < size; i++) {
    //     std::cout << Fxb_data[i] << " ";
    // }
    // std::cout << std::endl;

}

