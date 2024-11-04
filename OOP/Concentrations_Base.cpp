#include "Concentrations_Base.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"


Concentrations::Concentrations(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), EVol(mh.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi()), 
      gtPse(mesh_handler.GetTotalPse())

{
    
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();


}

void Concentrations::SetInitialValues(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation) {
    
    // Shared initialization tasks
    SetInitialConcentration(Cn, initial_value);

    if (perform_lithiation) {
        LithiationCalculation(Cn, psx);
    }

}

void Concentrations::SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value) {
    
    for (int i = 0; i < Cn.Size(); ++i) {
        Cn(i) = initial_value;
    }

}

void Concentrations::SetUpSolver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother) {
    
    M = new mfem::ParBilinearForm(fespace);

    Ps_gf = new mfem::ParGridFunction(fespace);
    *Ps_gf = psx;

    cP = new mfem::GridFunctionCoefficient(Ps_gf);

    M->AddDomainIntegrator(new mfem::MassIntegrator(*cP));
    M->Assemble();
    M->Finalize();

    mfem::HypreParMatrix HPM;
    M->FormSystemMatrix(boundary_dofs, HPM); // should be form linear system
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

void Concentrations::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {

    Rx2 = Rx1;
    Rx2 *= value;

}

void Concentrations::ForceTerm(mfem::ParGridFunction &gfc, mfem::ParLinearForm &Fxx, mfem::Array<int> boundary, mfem::ProductCoefficient m, bool apply_boundary_conditions) {

    std::unique_ptr<ParLinearForm> Bx2(new ParLinearForm(fespace));	

    Rxx = new ParGridFunction(fespace);
    *Rxx = gfc;

    cXx = new GridFunctionCoefficient(Rxx);

    Bx2->AddDomainIntegrator(new DomainLFIntegrator(*cXx));

    if (apply_boundary_conditions) {
        Bx2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m), boundary);
    }

    Bx2->Assemble();
    Fxx = std::move(*Bx2);

}

void Concentrations::TotalReaction(mfem::ParGridFunction &Rx, double value) {

    value = 0.0;
    Array<double> VtxVal(nC);
    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Rx.GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        value += EAvg(ei) * EVol(ei);
    }

    double global_value;
    MPI_Allreduce(&value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    infx = global_value / (mesh_handler.L_w);

}

std::shared_ptr<mfem::GridFunctionCoefficient> Concentrations::Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte ){

    mfem::ParGridFunction *Dx = new mfem::ParGridFunction(fespace);

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