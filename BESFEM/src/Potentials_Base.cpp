// /**
//  * @file Potentials_Base.cpp
//  * @brief Implementation of the potentials class methods for battery simulations.
//  */

#include "../include/Potentials_Base.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "mfem.hpp"

// /**
//  * @brief Constructor: Initializes the Potentials object.
//  */
PotentialBase::PotentialBase(Initialize_Geometry &geo, Domain_Parameters &para)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), EVol(para.EVol), px0(fespace.get()),
    TmpF(fespace.get())

{
    nE = geometry.nE; 
    nC = geometry.nC; 
    nV = geometry.nV; 

    X0 = mfem::HypreParVector(fespace.get()); // Initialize the potential vector
    px0 = mfem::ParGridFunction(fespace.get()); // Initialize the potential grid function
    TmpF = mfem::ParGridFunction(fespace.get()); // Temporary grid function for error calculations
}

// void PotentialBase::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);
// void PotentialBase::AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential);
// void PotentialBase::UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);


// void Potentials::SetInitialPotentials(mfem::ParGridFunction &ph, double initial_value) {
//     for (int i = 0; i < ph.Size(); ++i) {
//         ph(i) = initial_value; // Assign the initial value to each DoF
//     }
// }


// void Potentials::AssembleForceVector(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value, mfem::GridFunctionCoefficient &coef, std::unique_ptr<mfem::ParLinearForm> &rhs_form, mfem::ParLinearForm &rhs_form2)
// { 
//     Rx2 = Rx1;  // Copy the input field
//     Rx2 *= value; // Scale the field by the provided factor

//     coef.SetGridFunction(&Rx2); // Set the coefficient to the scaled field

//     rhs_form->Assemble();

//     rhs_form2 = *rhs_form; // Move the updated linear form to the right-hand side vector

// }



// void Potentials::ComputeGlobalError(mfem::ParGridFunction &px0, mfem::ParGridFunction &potential, mfem::ParGridFunction &psx, double &globalerror, double gtPsx)
// {
//     // Compute squared error using the auxiliary field
//     for (int vi = 0; vi < nV; vi++){
//         TmpF(vi) = pow(px0(vi) - potential(vi),2) * psx(vi);
//     }

//     double error_X = 0.0; // Initialize error accumulator
//     mfem::Array<double> VtxVal(nC);
//     mfem::Vector EAvg(nE);

//     // Calculate error contributions across all elements
//     for (int ei = 0; ei < nE; ei++){
//         TmpF.GetNodalValues(ei,VtxVal) ;
//         double val = 0.0;
//         for (int vt = 0; vt < nC; vt++){
//             val += VtxVal[vt];
//         }
//         EAvg(ei) = val/nC;	
//         error_X += EAvg(ei)*EVol(ei) ;					
//     }	
	
//     MPI_Allreduce(&error_X, &globalerror, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
//     globalerror /= gtPsx; // Normalize the error
//     globalerror = pow(globalerror, 0.5); // Compute the root mean square error

// }





// // void Potentials::SetUpSolver(mfem::CGSolver &solver, double value_1, double value_2) {
// //     solver.SetRelTol(value_1); // Set the relative tolerance
// //     solver.SetMaxIter(value_2); // Set the maximum number of iterations
// // }


// // void Potentials::PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix){
// //     smoother.SetType(mfem::HypreSmoother::Jacobi); // Set smoother type to Jacobi
// //     cg.SetPreconditioner(smoother); // Attach the preconditioner to the solver
// //     cg.SetOperator(KMatrix); // Set the stiffness matrix as the operator
// // }

// // void Potentials::ImplementBoundaryConditions(mfem::ConstantCoefficient &dbc_Coef, double Bv, mfem::ParGridFunction &phx, mfem::Array<int> dbc_bdr){
// //     dbc_Coef = mfem::ConstantCoefficient(Bv); // Set boundary value
// //     phx.ProjectBdrCoefficient(dbc_Coef, dbc_bdr); // Apply the boundary condition
// // }


// // void Potentials::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {
// //     Rx2 = Rx1;  // Copy the input field
// //     Rx2 *= value; // Scale the field by the provided factor
// // }


// // void Potentials::ForceVector(mfem::ParBilinearForm &K, mfem::Array<int> boundary, mfem::ParGridFunction &phx, 
// // mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B, mfem::ConstantCoefficient &Coef, mfem::Array<int> &bdr){

// //     phx.ProjectBdrCoefficient(Coef, bdr); // Apply boundary conditions
// //     K.FormLinearSystem(boundary, phx, plf_B, matrix, hpv_X, hpv_B); // Form the linear system
// // }


// // void Potentials::ErrorCalculation(mfem::ParGridFunction &phx, mfem::CGSolver &cg_solver, mfem::HypreParVector &fterm, mfem::ParGridFunction &psx, double error_X, double &globalerror, double gtPsx){

// //     *px0 = phx; // Store the current potential field
// //     px0->GetTrueDofs(*X0); // Extract degrees of freedom
// //     cg_solver.Mult(fterm, *X0); // Solve for the error term

// //     phx.Distribute(X0.get()); // Distribute the updated values
// //     mfem::ParGridFunction TmpF(fespace.get());

// //     // Compute squared error using the auxiliary field
// //     for (int vi = 0; vi < nV; vi++){
// //         TmpF(vi) = pow((*px0)(vi) - phx(vi),2) * psx(vi);
// //     }

// //     error_X = 0.0; // Initialize error accumulator
// //     mfem::Array<double> VtxVal(nC);
// //     mfem::Vector EAvg(nE);

// //     // Calculate error contributions across all elements
// //     for (int ei = 0; ei < nE; ei++){
// //         TmpF.GetNodalValues(ei,VtxVal) ;
// //         // double val = 0.0;
// //         // for (int vt = 0; vt < nC; vt++){
// //         //     val += VtxVal[vt];
// //         // }
// //         double val = std::accumulate(VtxVal.begin(), VtxVal.end(), 0.0);
// //         EAvg(ei) = val/nC;	
// //         error_X += EAvg(ei)*EVol(ei) ;					
// //     }	
	
// //     MPI_Allreduce(&error_X, &globalerror, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
// //     globalerror /= gtPsx; // Normalize the error
// //     globalerror = pow(globalerror, 0.5); // Compute the root mean square error

// // }