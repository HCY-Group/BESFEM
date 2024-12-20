/**
 * @file Concentrations_Base.cpp
 * @brief Implementation of the Concentrations class methods for battery simulations.
 */

#include "Concentrations_Base.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"

Concentrations::Concentrations(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), EVol(mh.EVol), gtPsi(mh.gtPsi), gtPse(mh.gtPse)

{
    
    nE = mesh_handler.nE; 
    nC = mesh_handler.nC; 
    nV = mesh_handler.nV; 


}

void Concentrations::SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value) {
    
    for (int i = 0; i < Cn.Size(); ++i) {
        Cn(i) = initial_value;
    }

}

void Concentrations::SetUpSolver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother) {
    
    // Create a parallel bilinear form for the finite element space.
    M = new mfem::ParBilinearForm(fespace);
    
    // Copy the provided potential field to a new grid function for internal operations.
    Ps_gf = new mfem::ParGridFunction(fespace);
    *Ps_gf = psx;

    // Wrap the grid function in a coefficient to use in the mass matrix integrator.
    cP = new mfem::GridFunctionCoefficient(Ps_gf);

    // Add a domain integrator to the bilinear form using the weighted mass integrator.
    M->AddDomainIntegrator(new mfem::MassIntegrator(*cP));
    
    // Assemble the bilinear form into a sparse matrix representation.
    M->Assemble();
    M->Finalize();

    // Construct the system matrix (mass matrix) and store it in HPM.
    mfem::HypreParMatrix HPM;
    M->FormSystemMatrix(boundary_dofs, HPM); // should be form linear system
    
    // Assign the constructed matrix to the shared pointer provided.
    Mmat = std::make_shared<mfem::HypreParMatrix>(HPM);

    // Configure the preconditioner using a Jacobi smoother.
    smoother.SetType(mfem::HypreSmoother::Jacobi);

    // Set up the solver for the mass matrix.
    m_solver.iterative_mode = false; // Use direct solving for the system matrix.
    m_solver.SetRelTol(1e-7); // Set relative tolerance for the solver.
    m_solver.SetAbsTol(0); // Set absolute tolerance for the solver.
    m_solver.SetMaxIter(102); // Limit the maximum number of iterations.
    m_solver.SetPrintLevel(0); // Suppress output from the solver.
    m_solver.SetPreconditioner(smoother); // Attach the preconditioner to the solver.
    m_solver.SetOperator(*Mmat); // Set the mass matrix as the operator to solve.

}

void Concentrations::LithiationCalculation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
    // Temporary grid function to store the product of concentration and potential.
    mfem::ParGridFunction TmpF(fespace);
    TmpF = Cn; // Copy concentration values to the temporary grid function.
    TmpF *= psx; // Element-wise multiply concentration by potential.


    double lSum = 0.0; // Local sum of lithiation contributions.
    mfem::Array<double> VtxVal(nC); // Array to store nodal values for each element.
    mfem::Vector EAvg(nE); // Vector to store the average contribution for each element.

    for (int ei = 0; ei < nE; ei++) {
        TmpF.GetNodalValues(ei, VtxVal); // Extract nodal values for the current element.
        double val = 0.0; // Temporary variable to accumulate nodal contributions.
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC; // Compute the average contribution for the element.     
        lSum += EAvg(ei) * EVol(ei); // Weight by element volume and add to the local sum.
    }

    double gSum; // Global sum to aggregate contributions across all MPI processes.
    MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Perform global reduction.
    double Xfr = gSum / gtPsi; // Calculate the degree of lithiation as the normalized sum.
}

void Concentrations::ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF) {
    PGF = psx; // Copy the input potential field to the target grid function.
    PGF.Neg(); // Negate all values in the target grid function.
}

void Concentrations::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {

    Rx2 = Rx1; // Copy the input reaction field to the output reaction field.
    Rx2 *= value; // Scale the output reaction field by the specified factor.

}


void Concentrations::ForceTerm(mfem::ParGridFunction &gfc, mfem::ParLinearForm &Fxx, mfem::Array<int> boundary, mfem::ProductCoefficient m, bool apply_boundary_conditions) {
    
    // Create a unique pointer to a parallel linear form associated with the finite element space.
    std::unique_ptr<mfem::ParLinearForm> Bx2(new mfem::ParLinearForm(fespace));	
    
    // Initialize a new ParGridFunction to hold the input field and copy values from gfc.
    Rxx = new mfem::ParGridFunction(fespace);
    *Rxx = gfc;
    
    // Create a GridFunctionCoefficient from the ParGridFunction for use in integrators.
    cXx = new mfem::GridFunctionCoefficient(Rxx);

    // Add a domain integrator to compute contributions from the entire domain.
    Bx2->AddDomainIntegrator(new DomainLFIntegrator(*cXx));

    // If boundary conditions are to be applied, add a boundary integrator.
    if (apply_boundary_conditions) {
        Bx2->AddBoundaryIntegrator(new BoundaryLFIntegrator(m), boundary);
    }
    
    // Assemble the linear form to finalize its representation.
    Bx2->Assemble();
    
    // Move the assembled linear form into the provided Fxx reference.
    Fxx = std::move(*Bx2);

}


void Concentrations::TotalReaction(mfem::ParGridFunction &Rx, double xCrnt) {
    
    xCrnt = 0.0; // Initialize the local total reaction value to zero.
    mfem::Array<double> VtxVal(nC); // Array to hold nodal values for each element.
    mfem::Vector EAvg(nE); // Vector to store average reaction contributions for each element.

    for (int ei = 0; ei < nE; ei++) {
        Rx.GetNodalValues(ei, VtxVal); // Retrieve the nodal values of the reaction field for the current element.
        double val = 0.0; // Temporary variable to accumulate nodal contributions.
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC; // Calculate the average reaction value for the element.
        xCrnt += EAvg(ei) * EVol(ei); // Weight by the element volume and add to the total reaction.
    }
    
    // Perform a global reduction to sum up contributions across all MPI processes.
    MPI_Allreduce(&xCrnt, &geCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // Calculate the reaction current density by normalizing with the domain characteristic length.
    infx = geCrnt / (mesh_handler.L_w);

}


std::shared_ptr<mfem::GridFunctionCoefficient> Concentrations::Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte ){
    
    // Create a new parallel grid function to store the computed diffusivity values.
    mfem::ParGridFunction *Dx = new mfem::ParGridFunction(fespace);
    
    // Loop through all vertices in the domain to calculate diffusivity.
    for (int vi = 0; vi < nV; vi++) {
        if (particle_electrolyte) {
            
            // Compute diffusivity for the particle based on a polynomial model.
            (*Dx)(vi) = psx(vi) * (0.0277 - 0.084 * Cn(vi) + 0.1003 * Cn(vi) * Cn(vi)) * 1.0e-8;
            
            // Cap the diffusivity at a maximum value for stability.
            if ((*Dx)(vi) > 4.6e-10) {
                (*Dx)(vi) = 4.6e-10;
            }

        } else {
            // Compute diffusivity for the electrolyte based on an exponential model.
            (*Dx)(vi) = psx(vi) * Constants::D0 * exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
        }
    }
    
    // Wrap the computed diffusivity grid function in a GridFunctionCoefficient and return it.
    return std::make_shared<mfem::GridFunctionCoefficient>(Dx);

}


void Concentrations::KMatrix(mfem::Array<int> boundary, mfem::ParGridFunction &Cn, mfem::ParLinearForm &Fxx, std::shared_ptr<mfem::HypreParMatrix> &Kmatx, mfem::HypreParVector &X1v, mfem::HypreParVector &Fxb, mfem::GridFunctionCoefficient *cDx) {
    
    // Create a new unique pointer to a parallel bilinear form for the finite element space.
    std::unique_ptr<mfem::ParBilinearForm> Kx2(new mfem::ParBilinearForm(fespace));
    
    // Temporary matrix to hold the assembled stiffness matrix.
    mfem::HypreParMatrix Khpm;
    
    // Add a domain integrator for the diffusion term using the given diffusivity coefficient (cDx).
    Kx2->AddDomainIntegrator(new DiffusionIntegrator(*cDx));
    
    // Assemble the bilinear form into a sparse matrix.
    Kx2->Assemble();
    
    // Form the linear system, considering boundary conditions, concentration field (Cn), and the right-hand side (Fxx).
    Kx2->FormLinearSystem(boundary, Cn, Fxx, Khpm, X1v, Fxb);
    
    // Convert the assembled stiffness matrix into a shared pointer for further use.
    Kmatx = std::make_shared<mfem::HypreParMatrix>(Khpm);
    
    // Scale the right-hand side vector for time-stepping by multiplying with the time step constant (Constants::dt).
    Fxb *= Constants::dt;

}


void Concentrations::SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {

    CeC = 0.0; // Initialize total salt concentration to zero.
    
    // Temporary grid function to hold the product of concentration and potential fields.
    mfem::ParGridFunction CeT(fespace);
    mfem::Array<double> VtxVal(nC); // Array to store nodal values for each element.
    mfem::Vector EAvg(nE); // Vector to store average concentration values for each element.

    // Calculate the product of concentration and potential for each element.
    CeT = Cn;
    CeT *= psx;
    
    // Loop through all elements to calculate the average concentration for each element.
    for (int ei = 0; ei < nE; ei++){
        CeT.GetNodalValues(ei,VtxVal); // Retrieve the nodal values for the current element. 
        double val = 0.0; // Temporary variable to accumulate nodal contributions.
        
        for (int vt = 0; vt < nC; vt++){
            val += VtxVal[vt];
        }

        EAvg(ei) = val/nC; // Compute the average concentration for the element.
        
        // Update the total salt concentration by adding the weighted average for the element.
        CeC += EAvg(ei)*EVol(ei); 
    }
    
    // Perform a global reduction to sum up the salt concentration across all MPI processes.
    MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
    
    // Compute the average concentration across the entire electrolyte.
    CeAvg = gCeC/gtPse;	
    
    // Adjust the concentration field by normalizing the average concentration with the initial value.
    Cn -= (CeAvg-Ce0);
    MPI_Barrier(MPI_COMM_WORLD);

}

