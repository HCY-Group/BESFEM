#include "CnP.hpp"
#include "mfem.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

void InitializeCnP(mfem::ParFiniteElementSpace &fespace,
                mfem::ParGridFunction &psi,
                double gtPsi,
                mfem::ParMesh &pmesh)
{
    // initial condition
    ParGridFunction CnP(&fespace);
    double Cp0 = 0.3; // initial value
    CnP = Cp0;

    // degree of lithiation
    double Xfr = 0.0;

    ParGridFunction TmpF(&fespace);
    TmpF = CnP;
    TmpF *= psi;
    double lSum = 0.0;
    int nE = pmesh.GetNE();
    int nC = pow(2, pmesh.Dimension());
    Vector EVol(nE);
    for (int ei = 0; ei < nE; ei++) {
        EVol(ei) = pmesh.GetElementVolume(ei);
    }

    Vector EAvg(nE);
    Vector VtxVal(nC);
    Array<double> nodalValues;
    for (int ei = 0; ei < nE; ei++) {
        TmpF.GetNodalValues(ei, nodalValues);
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += nodalValues[vt];
        }
        EAvg(ei) = val / nC;
        lSum += EAvg(ei) * EVol(ei);
    }

    cout << "Sum: " << lSum << std::endl;


    double gSum;
    MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Xfr = gSum / gtPsi;

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
    HypreParMatrix Mmatp;
    std::unique_ptr<ParBilinearForm> Mt(new ParBilinearForm(&fespace));
    GridFunctionCoefficient cPs(&psi);
    Mt->AddDomainIntegrator(new MassIntegrator(cPs));
    Mt->Assemble();
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

    HypreParMatrix *Tmatp;

    // SBM stiffness matrix
    ParGridFunction Dp(&fespace);

    HypreParMatrix Kmatp;

    // force vector
    ParGridFunction Rxc(&fespace);
    ParLinearForm *Bc2;
    ParLinearForm Fct(&fespace);
    HypreParVector Fcb(&fespace);

    // create a Vector for CnP
    HypreParVector CpV0(&fespace), CpVn(&fespace), RHCp(&fespace);

    int nDof = CpV0.Size();

    // Vector of psi
    HypreParVector PsVc(&fespace);
    psi.GetTrueDofs(PsVc);

    cout << "Degree of Lithiation: " << Xfr << std::endl;
}
