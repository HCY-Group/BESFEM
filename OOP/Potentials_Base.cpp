/**
 * @file Potentials_Base.cpp
 * @brief Implementation of the potentials class methods for battery simulations.
 */

#include "Potentials_Base.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"

/**
 * @brief Constructor: Initializes the Potentials object.
 */
Potentials::Potentials(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), EVol(mh.EVol)

{
    
    nE = mesh_handler.nE; 
    nC = mesh_handler.nC; 
    nV = mesh_handler.nV; 

    // px0 = new mfem::ParGridFunction(fespace); // Potential values before iteration
    px0 = std::make_unique<mfem::ParGridFunction>(fespace);
    Rxx = std::make_unique<mfem::ParGridFunction>(fespace);
    // Bx2 = std::make_unique<mfem::ParLinearForm>(fespace);

    X0 = mfem::HypreParVector(fespace); // Hypre vector for solving linear systems

}


void Potentials::SetInitialPotentials(mfem::ParGridFunction &ph, double initial_value) {
    
    for (int i = 0; i < ph.Size(); ++i) {
        ph(i) = initial_value; // Assign the initial value to each DoF
    }

}

void Potentials::SetUpSolver(mfem::CGSolver &solver, double value_1, double value_2) {
    
    solver.SetRelTol(value_1); // Set the relative tolerance
    solver.SetMaxIter(value_2); // Set the maximum number of iterations

}


void Potentials::KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, 
mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B){

    K.Update(); // Update the bilinear form
    K.AddDomainIntegrator(new DiffusionIntegrator(gfc)); // Add domain integrator for stiffness
    K.Assemble(); // Assemble the matrix
    K.FormLinearSystem(boundary, potential, plf_B, matrix, hpv_X, hpv_B); // Form the linear system
}

void Potentials::PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix){

    smoother.SetType(HypreSmoother::Jacobi); // Set smoother type to Jacobi
    cg.SetPreconditioner(smoother); // Attach the preconditioner to the solver
    cg.SetOperator(KMatrix); // Set the stiffness matrix as the operator


}

void Potentials::ImplementBoundaryConditions(mfem::ConstantCoefficient &dbc_Coef, double Bv, mfem::ParGridFunction &phx, mfem::Array<int> dbc_bdr){

    dbc_Coef = mfem::ConstantCoefficient(Bv); // Set boundary value
    phx.ProjectBdrCoefficient(dbc_Coef, dbc_bdr); // Apply the boundary condition


}


void Potentials::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {

    Rx2 = Rx1;  // Copy the input field
    Rx2 *= value; // Scale the field by the provided factor

}

void Potentials::ForceTerm(mfem::ParGridFunction &Rx2, mfem::ParLinearForm &Fxx) {

    Bx2 = std::make_unique<mfem::ParLinearForm>(fespace);
    // Bx2->Update();

    // Rxx = new mfem::ParGridFunction(fespace); // Initialize intermediate reaction field
    *Rxx = Rx2;

    // cXx = new mfem::GridFunctionCoefficient(Rxx.get()); // Wrap as a coefficient
    cXx = std::make_unique<mfem::GridFunctionCoefficient>(Rxx.get());

    Bx2->AddDomainIntegrator(new mfem::DomainLFIntegrator(*cXx)); // Integrate over the domain

    Bx2->Assemble(); // Assemble the linear form
    Fxx = std::move(*Bx2); // Transfer ownership to the provided reference

}


void Potentials::ForceVector(mfem::ParBilinearForm &K, mfem::Array<int> boundary, mfem::ParGridFunction &phx, 
mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B, mfem::ConstantCoefficient &Coef, mfem::Array<int> &bdr){

    phx.ProjectBdrCoefficient(Coef, bdr); // Apply boundary conditions
    K.FormLinearSystem(boundary, phx, plf_B, matrix, hpv_X, hpv_B); // Form the linear system
}



void Potentials::ErrorCalculation(mfem::ParGridFunction &phx, mfem::CGSolver &cg_solver, mfem::HypreParVector &fterm, mfem::ParGridFunction &psx, double error_X, double &globalerror, double gtPsx){

    *px0 = phx; // Store the current potential field
    px0->GetTrueDofs(X0); // Extract degrees of freedom
    cg_solver.Mult(fterm, X0); // Solve for the error term

    phx.Distribute(X0); // Distribute the updated values
    mfem::ParGridFunction TmpF(fespace);

    // Compute squared error using the auxiliary field
    for (int vi = 0; vi < nV; vi++){
        TmpF(vi) = pow((*px0)(vi) - phx(vi),2) * psx(vi);
    }

    error_X = 0.0; // Initialize error accumulator
    mfem::Array<double> VtxVal(nC);
    mfem::Vector EAvg(nE);

    // Calculate error contributions across all elements
    for (int ei = 0; ei < nE; ei++){
        TmpF.GetNodalValues(ei,VtxVal) ;
        // double val = 0.0;
        // for (int vt = 0; vt < nC; vt++){
        //     val += VtxVal[vt];
        // }
        double val = std::accumulate(VtxVal.begin(), VtxVal.end(), 0.0);
        EAvg(ei) = val/nC;	
        error_X += EAvg(ei)*EVol(ei) ;					
    }	
	
    MPI_Allreduce(&error_X, &globalerror, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
    globalerror /= gtPsx; // Normalize the error
    globalerror = pow(globalerror, 0.5); // Compute the root mean square error

}