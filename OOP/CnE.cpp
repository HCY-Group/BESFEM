#include "CnE.hpp"
#include "mfem.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

void InitializeCnE(mfem::ParFiniteElementSpace &fespace,
                mfem::ParGridFunction &pse,
                mfem::ParMesh &pmesh)
{
    // initial condition
    ParGridFunction CnE(&fespace);
    double Ce0 = 0.001; // initial value
    CnE = Ce0;

    // Some containers that will be used later.
    Array<int> boundary_dofs; // nature boundary

    // Neumann BC on the west boundary. CnE
    Array<int> nbc_w_bdr(pmesh.bdr_attributes.Max());
    nbc_w_bdr = 0;
    nbc_w_bdr[0] = 1;

    // Dirichlet BC on the east boundary. phP
    Array<int> dbc_e_bdr(pmesh.bdr_attributes.Max());
    dbc_e_bdr = 0;
    dbc_e_bdr[2] = 1;
    // use dbc_e_bdr array to extract all node labels of Dirichlet BC
    Array<int> ess_tdof_list_e(0);
    fespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

    // Dirichlet BC on the west boundary. phE
    Array<int> dbc_w_bdr(pmesh.bdr_attributes.Max());
    dbc_w_bdr = 0;
    dbc_w_bdr[0] = 1;
    // use dbc_w_bdr array to extract all node labels of Dirichlet BC
    Array<int> ess_tdof_list_w(0);
    fespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

    // SBM mass matrix
    HypreParMatrix Mmate;
    std::unique_ptr<ParBilinearForm> Me(new ParBilinearForm(&fespace));
    GridFunctionCoefficient cPe(&pse);
    Me->AddDomainIntegrator(new MassIntegrator(cPe));
    Me->Assemble();
    Me->FormSystemMatrix(boundary_dofs, Mmate);

    HypreSmoother Me_prec;
    CGSolver Me_solver(MPI_COMM_WORLD);

    Me_solver.iterative_mode = false;
    Me_solver.SetRelTol(1e-7);
    Me_solver.SetAbsTol(0);
    Me_solver.SetMaxIter(200);
    Me_solver.SetPrintLevel(0);
    Me_prec.SetType(HypreSmoother::Jacobi);
    Me_solver.SetPreconditioner(Me_prec);
    Me_solver.SetOperator(Mmate);

    HypreParMatrix *TmatL, *TmatR;

    // stiffness matrix
    ParGridFunction De(&fespace);
    HypreParMatrix Kmate;
    // std::unique_ptr<ParBilinearForm> Ke2(new ParBilinearForm(&fespace));
    // Ke2->AddDomainIntegrator(new DiffusionIntegrator(cPe));
    // Ke2->Assemble();
    // Ke2->FormSystemMatrix(boundary_dofs, Kmate);

    // force vector
    ParGridFunction Rxe(&fespace);
    // std::unique_ptr<ParLinearForm> Be2(new ParLinearForm(&fespace));
    // Be2->AddDomainIntegrator(new DomainLFIntegrator(cPe));
    // Be2->Assemble();
    // ParLinearForm Fet(&fespace);
    // Fet = std::move(*Be2);

    ParLinearForm *Be2;
    ParLinearForm Fet(&fespace);
    HypreParVector Feb(&fespace);


    // Neumann BC on the west boundary
    ParGridFunction PeR(&fespace);
    PeR = pse;
    PeR.Neg();
    GridFunctionCoefficient matCoef_R(&PeR);

    // Vectors for CnE
    HypreParVector CeV0(&fespace), CeVn(&fespace), RHSe(&fespace);

    // parameters used in the calculations
    double eCrnt = 0.0;
    double geCrnt = 0.0;
    double infx = 0.0;

    ParGridFunction CeT(&fespace);
    double CeC = 0.0;
    double CeAvg = 0.0;
    double gCeC = 0.0;

    // Compute average value
    // ParGridFunction TmpF(&fespace);
    // TmpF = CnE;
    // TmpF *= pse;
    // double lSum = 0.0;
    // int nE = pmesh.GetNE();
    // int nC = pow(2, pmesh.Dimension());
    // Vector EVol(nE);
    // for (int ei = 0; ei < nE; ei++) {
    //     EVol(ei) = pmesh.GetElementVolume(ei);
    // }

    // Vector EAvg(nE);
    // Vector VtxVal(nC);
    // Array<double> nodalValues;
    // for (int ei = 0; ei < nE; ei++) {
    //     TmpF.GetNodalValues(ei, nodalValues);
    //     double val = 0.0;
    //     for (int vt = 0; vt < nC; vt++) {
    //         val += nodalValues[vt];
    //     }
    //     EAvg(ei) = val / nC;
    //     lSum += EAvg(ei) * EVol(ei);
    // }

    // double gSum;
    // MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // CeAvg = gSum / gtPsi;

    // cout << "Average Concentration of Species E: " << CeAvg << std::endl;
}
