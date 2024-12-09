#include "Potentials_Base.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"


Potentials::Potentials(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), EVol(mh.GetElementVolume())

{
    
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();

    px0 = new mfem::ParGridFunction(fespace); // values before iteration
    X0 = mfem::HypreParVector(fespace);

    // TmpF = mfem::ParGridFunction(fespace);




}


void Potentials::SetInitialPotentials(mfem::ParGridFunction &ph, double initial_value) {
    
    for (int i = 0; i < ph.Size(); ++i) {
        ph(i) = initial_value;
    }

}

void Potentials::SetUpSolver(mfem::CGSolver &solver, double value_1, double value_2) {
    
    solver.SetRelTol(value_1);
    solver.SetMaxIter(value_2);

}


void Potentials::KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, 
mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B){

    K.Update();
    K.AddDomainIntegrator(new DiffusionIntegrator(gfc));
    K.Assemble();
    K.FormLinearSystem(boundary, potential, plf_B, matrix, hpv_X, hpv_B);
}

void Potentials::PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix){

    smoother.SetType(HypreSmoother::Jacobi);
    cg.SetPreconditioner(smoother);
    cg.SetOperator(KMatrix);


}

void Potentials::ImplementBoundaryConditions(mfem::ConstantCoefficient &dbc_Coef, double Bv, mfem::ParGridFunction &phx, mfem::Array<int> dbc_bdr){

    dbc_Coef = mfem::ConstantCoefficient(Bv);
    phx.ProjectBdrCoefficient(dbc_Coef, dbc_bdr);


}


void Potentials::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {

    Rx2 = Rx1;
    Rx2 *= value;

}

void Potentials::ForceTerm(mfem::ParGridFunction &Rx2, mfem::ParLinearForm &Fxx) {

    std::unique_ptr<ParLinearForm> Bx2(new ParLinearForm(fespace));	

    Rxx = new mfem::ParGridFunction(fespace);
    *Rxx = Rx2;
    cXx = new mfem::GridFunctionCoefficient(Rxx);

    Bx2->AddDomainIntegrator(new DomainLFIntegrator(*cXx));
    Bx2->Assemble();
    Fxx = std::move(*Bx2);

}


void Potentials::ForceVector(mfem::ParBilinearForm &K, mfem::Array<int> boundary, mfem::ParGridFunction &phx, 
mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B, mfem::ConstantCoefficient &Coef, mfem::Array<int> &bdr){

    phx.ProjectBdrCoefficient(Coef, bdr);
    K.FormLinearSystem(boundary, phx, plf_B, matrix, hpv_X, hpv_B);

}


void Potentials::ErrorCalculation(mfem::ParGridFunction &phx, mfem::CGSolver &cg_solver, mfem::HypreParVector &fterm, mfem::ParGridFunction &psx, double error_X, double &globalerror, double gtPsx){

    *px0 = phx;
    px0->GetTrueDofs(X0);
    cg_solver.Mult(fterm, X0);

    phx.Distribute(X0);
    mfem::ParGridFunction TmpF(fespace);


    for (int vi = 0; vi < nV; vi++){
        TmpF(vi) = pow((*px0)(vi) - phx(vi),2) * psx(vi);
    }

    error_X = 0.0;
    mfem::Array<double> VtxVal(nC);
    mfem::Vector EAvg(nE);


    for (int ei = 0; ei < nE; ei++){
        TmpF.GetNodalValues(ei,VtxVal) ;
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++){
            val += VtxVal[vt];
        }
        EAvg(ei) = val/nC;	
        error_X += EAvg(ei)*EVol(ei) ;					
    }	
	
    MPI_Allreduce(&error_X, &globalerror, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
    globalerror /= gtPsx;
    globalerror = pow(globalerror, 0.5);

    // std::cout << "global error: " << globalerror_X << std::endl;
    // std::cout << " error: " << error_X << std::endl;



}