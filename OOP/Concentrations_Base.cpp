#include "Concentrations_Base.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"


Concentrations::Concentrations(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), Mmat(std::make_shared<mfem::HypreParMatrix>()), solver(std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD)),
    EVol(mh.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi()), gtPse(mesh_handler.GetTotalPse())

{
    
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();


}

void Concentrations::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother, bool perform_lithiation) {
    
    // Shared initialization tasks
    CreateCn(Cn, initial_value, perform_lithiation, psx);
    Solver(psx, Mmat, m_solver, smoother);

    // Class-specific tasks defined in inherit classes
    Initialize_Differences(psx);
}

void Concentrations::CreateCn(mfem::ParGridFunction &Cn, double initial_value, bool perform_lithiation, mfem::ParGridFunction &psx) {
    for (int i = 0; i < Cn.Size(); ++i) {
        Cn(i) = initial_value;
    }
    if (perform_lithiation) {
        LithiationCalculation(Cn, psx);
    }
}

void Concentrations::Solver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother) {
    
    M = new mfem::ParBilinearForm(fespace);

    Ps_gf = new mfem::ParGridFunction(fespace);
    *Ps_gf = psx;

    cP = new mfem::GridFunctionCoefficient(Ps_gf);

    M->AddDomainIntegrator(new mfem::MassIntegrator(*cP));
    M->Assemble();
    M->Finalize();

    mfem::HypreParMatrix HPM;
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

void Concentrations::LithiationCalculation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
    mfem::ParGridFunction TmpF(fespace);
    TmpF = Cn;
    TmpF *= psx;


    double lSum = 0.0;
    mfem::Array<double> VtxVal(nC);
    mfem::Vector EAvg(nE);
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

    // std::cout << "gSum: " << gSum << std::endl;


}

void Concentrations::ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF) {
    PGF = psx;
    PGF.Neg();
}
