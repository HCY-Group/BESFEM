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

Concentrations::Concentrations(MeshHandler &mesh_handler)
    : mesh_handler(mesh_handler), fespace(mesh_handler.GetFESpace()), pmesh(mesh_handler.GetPmesh()), psi(*mesh_handler.GetPsi()), pse(*mesh_handler.GetPse()),
      EVol(mesh_handler.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi()), reaction(mesh_handler, *this), PeR(fespace.get()), matCoef_R(&PeR)

      
{

    nbc_w_bdr.SetSize(mesh_handler.pmesh->bdr_attributes.Max());
    nbc_w_bdr = 0;  // Initialize the array with zeros

    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();
    
    // Initialize ParGridFunction for CnP and CnE
    CnP = make_unique<ParGridFunction>(fespace.get());
    CnE = make_unique<ParGridFunction>(fespace.get());

    Rxc = make_unique<ParGridFunction>(fespace.get());
    Rxe = make_unique<ParGridFunction>(fespace.get());

    X1v = HypreParVector(fespace.get());
    Fcb = HypreParVector(fespace.get());
    Feb = HypreParVector(fespace.get());

    // Mmatp = std::make_shared<mfem::HypreParMatrix>();
    // Mmate = std::make_shared<mfem::HypreParMatrix>();

    reaction.Initialize(); 
}


void Concentrations::InitializeCnP(std::shared_ptr<ParFiniteElementSpace> fespace) {

    Lithiation(*CnP, 0.3, fespace);

    Mmatp = std::make_shared<mfem::HypreParMatrix>();

    SBM_Matrix(psi, Mmatp, fespace);

    Solver(Mmatp);

    // CnP->Print(std::cout);

}

void Concentrations::InitializeCnE(std::shared_ptr<ParFiniteElementSpace> fespace) {

    CreateCnE(*CnE, 0.001);

    // Mmate = new HypreParMatrix();
    // SetupBoundaryConditions();
    Mmate = std::make_shared<HypreParMatrix>();

    SBM_Matrix(pse, Mmate, fespace);

    // std::cout << "Mmate Rows after SBM in Initialize: " << Mmate->GetGlobalNumRows() << std::endl;
    // std::cout << "Mmate Columns after SBM in Initialize: " << Mmate->GetGlobalNumCols() << std::endl;

    Solver(Mmate);
    ImposeNeumannBC(PeR, pse);
    // GridFunctionCoefficient matCoef_R(&PeR);

    // PeR.Print(std::cout);

}

void Concentrations::TimeStepCnP(std::shared_ptr<ParFiniteElementSpace> fespace) {

    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    GridFunctionCoefficient cAp(Rxc.get());
    SetupRx(Rxn, *Rxc, Constants::rho, cAp); // Rxn needs to come from Reacions.cpp

    Array<int> dummy_boundary;
    
    ConstantCoefficient dummy_coef(0.0);

    ForceTerm(fespace, cAp, Fct, dummy_boundary, dummy_coef, false); // false since not applying BCs

    // GridFunctionCoefficient cDp = Diffusivity(psi, *CnP, true);
    std::shared_ptr<GridFunctionCoefficient> cDp = Diffusivity(psi, *CnP, true); // true since using first equation

    // std::cout << "fespace address in TimeStepCnP: " << fespace.get() << std::endl;
    mfem::HypreParMatrix Kmatp;
    K_Matrix(boundary_dofs, *CnP, Fct, Kmatp, X1v, Fcb, cDp.get());

    // std::cout << "Mmatp Rows in Timestep: " << Mmatp->GetGlobalNumRows() << std::endl;
    // // std::cout << "Kmatp Rows: " << Kmatp.GetGlobalNumRows() << std::endl;
    // std::cout << "Mmatp Columns in Timestep: " << Mmatp->GetGlobalNumCols() << std::endl;
    // // std::cout << "Kmatp Columns: " << Kmatp.GetGlobalNumCols() << std::endl;
    
    
    // assert(Mmatp.GetGlobalNumRows() == Kmatp.GetGlobalNumRows());
    // assert(Mmatp.GetGlobalNumCols() == Kmatp.GetGlobalNumCols());
    // Create T Matrix
    // Tmatp = 1.0 * Mmatp + dt * Kmatp
    // Tmatp = Add(1.0, Mmatp, Constants::dt, Kmatp); // should be -dt;

    // Degree of Lithiation
    // LithiationCalculation(*CnP, fespace);

}

void Concentrations::TimeStepCnE(std::shared_ptr<ParFiniteElementSpace> fespace) {

    mfem::ParGridFunction &Rxn = *reaction.Rxn;
    GridFunctionCoefficient cAe(Rxe.get());
    SetupRx(Rxn, *Rxe, Constants::t_minus, cAe);

    TotalReaction(*Rxe, eCrnt);
    ConstantCoefficient nbcCoef(-infx); // Neumann BC works with -infx

    ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    ForceTerm(fespace, cAe, Fet, nbc_w_bdr, nbcCoef, true); // true since applying boundary conditions

    std::shared_ptr<GridFunctionCoefficient> cDe = Diffusivity(pse, *CnE, false); // false using other equation

    // std::cout << "fespace address in TimeStepCnE: " << fespace.get() << std::endl;

    HypreParMatrix Kmate;
    K_Matrix(boundary_dofs, *CnE, Fet, Kmate, X1v, Feb, cDe.get());

    // std::cout << "Mmate Rows in Timestep: " << Mmate->GetGlobalNumRows() << std::endl;
    // // std::cout << "Kmatp Rows: " << Kmatp.GetGlobalNumRows() << std::endl;
    // std::cout << "Mmate Columns in Timestep: " << Mmate->GetGlobalNumCols() << std::endl;

}

void Concentrations::CreateCnE(mfem::ParGridFunction &Cn, double initial_value) {

    Cn = initial_value;

}

void Concentrations::LithiationCalculation(mfem::ParGridFunction &Cn, std::shared_ptr<ParFiniteElementSpace> fespace) {

    // Cn = initial_value;

    ParGridFunction TmpF(fespace.get());
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

void Concentrations::Lithiation(mfem::ParGridFunction &Cn, double initial_value, std::shared_ptr<ParFiniteElementSpace> fespace) {

    Cn = initial_value;
    LithiationCalculation(Cn, fespace);
}

void Concentrations::SBM_Matrix(mfem::ParGridFunction &psx, std::shared_ptr<HypreParMatrix> &Mmat, std::shared_ptr<ParFiniteElementSpace> fespace) {

    // SetupBoundaryConditions();
    std::cout << "fespace address in SBM_Matrix: " << fespace.get() << std::endl;


    std::unique_ptr<ParBilinearForm> M(new ParBilinearForm(fespace.get()));
    GridFunctionCoefficient cP(&psx);
    M->AddDomainIntegrator(new MassIntegrator(cP));
    M->Assemble();
    M->Finalize();

    HypreParMatrix HPM;
    M->FormSystemMatrix(boundary_dofs, HPM);
    Mmat = std::make_shared<mfem::HypreParMatrix>(HPM);

    // std::cout << "Mmat Rows after assembly: " << Mmat->GetGlobalNumRows() << std::endl;
    // std::cout << "Mmat Columns after assembly: " << Mmat->GetGlobalNumCols() << std::endl;

}

// want to use the Mmatp here as was used in SBM function
void Concentrations::Solver(std::shared_ptr<HypreParMatrix> &Mmat) {
    
    std::cout << "fespace address in Solver: " << fespace.get() << std::endl;

    HypreSmoother M_prec;
    CGSolver M_solver(MPI_COMM_WORLD);

    M_solver.iterative_mode = false;
    M_solver.SetRelTol(1e-7);
    M_solver.SetAbsTol(0);
    M_solver.SetMaxIter(102);
    M_solver.SetPrintLevel(0);
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(*Mmat);


}

void Concentrations::SetupBoundaryConditions(std::shared_ptr<ParFiniteElementSpace> fespace) {
    
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

void Concentrations::ForceTerm(std::shared_ptr<ParFiniteElementSpace> fespace, GridFunctionCoefficient cXx, mfem::ParLinearForm &Fxx, Array<int> boundary, ConstantCoefficient m, bool apply_boundary_conditions) {

    std::unique_ptr<ParLinearForm> Bx2(new ParLinearForm(fespace.get()));	

    Bx2->AddDomainIntegrator(new DomainLFIntegrator(cXx));

    if (apply_boundary_conditions) {
        Bx2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m), boundary);
    }

    Bx2->Assemble();
    Fxx = std::move(*Bx2);

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

    mfem::ParGridFunction* Dx = new mfem::ParGridFunction(fespace.get());

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

void Concentrations::K_Matrix(Array<int> boundary, mfem::ParGridFunction &Cn, ParLinearForm &Fxx, HypreParMatrix &Kmatx, HypreParVector &X1v, HypreParVector &Fxb, GridFunctionCoefficient *cDx) {


    // std::cout << "Creating ParBilinearForm for fespace with size: " << fespace->GetTrueVSize() << std::endl;
    // std::cout << "Boundary DOFs size: " << boundary.Size() << std::endl;

    // SetupBoundaryConditions();
    std::unique_ptr<ParBilinearForm> Kx2(new ParBilinearForm(fespace.get()));
    // std::cout << "fespace address in K_Matrix: " << fespace.get() << std::endl;


    // std::cout << "cDx values in K_Matrix:" << std::endl;
    // for (int vi = 0; vi < fespace->GetTrueVSize(); ++vi) {
    //     std::cout << (*cDx->GetGridFunction())(vi) << " ";
    // }
    // std::cout << std::endl;
    
    Kx2->AddDomainIntegrator(new DiffusionIntegrator(*cDx));
    Kx2->Assemble();
    Kx2->FormLinearSystem(boundary, Cn, Fxx, Kmatx, X1v, Fxb);

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

