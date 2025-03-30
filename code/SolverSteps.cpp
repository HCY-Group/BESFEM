#include "SolverSteps.hpp"
#include "Initialize_Geometry.hpp"
#include "Constants.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

// // Constructor
// SolverSteps::SolverSteps(std::shared_ptr<mfem::ParFiniteElementSpace> fespace)
// : fespace(fespace)

// {}

// Constructor for shared_ptr
SolverSteps::SolverSteps(std::shared_ptr<mfem::ParFiniteElementSpace> fespace)
: fespace(fespace), raw_fespace(nullptr) {}

// Constructor for raw pointer
SolverSteps::SolverSteps(mfem::ParFiniteElementSpace* fespace)
: raw_fespace(fespace), fespace(nullptr) {}

// Destructor
SolverSteps::~SolverSteps() {}

// void SolverSteps::MassMatrix(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat) {
//     if (!fespace) {
//         throw std::runtime_error("fespace is null");
//     }
//     // Create a parallel bilinear form for the finite element space
//     M = new mfem::ParBilinearForm(fespace.get());

//     // Copy the provided potential field to a new grid function for internal operations
//     temp_ps = new mfem::ParGridFunction(fespace.get());
//     *temp_ps = psx;

//     // Wrap the grid function in a coefficient (psi value) to use in the mass matrix integrator
//     coef = new mfem::GridFunctionCoefficient(temp_ps);

//     // Add a domain integrator to the bilinear form using the weighted mass integrator
//     M->AddDomainIntegrator(new mfem::MassIntegrator(*coef));
    
//     // Assemble the bilinear form into a sparse matrix representation
//     M->Assemble();

//     // Construct the mass matrix and store it in HPM
//     mfem::HypreParMatrix HPM;
//     M->FormSystemMatrix(boundary_dofs, HPM); // should be form linear system
    
//     // Assign the mass matrix to the shared pointer
//     Mmat = std::make_shared<mfem::HypreParMatrix>(HPM);

//     // // Clean up dynamically allocated objects
//     // delete coef; 
//     // delete temp_ps;
//     // delete M;

// };

// void SolverSteps::StiffnessMatrix(std::shared_ptr<mfem::GridFunctionCoefficient> cDx, mfem::Array<int> boundary, mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, std::shared_ptr<mfem::HypreParMatrix> &Kmat, mfem::HypreParVector &RHS){

//     // Create a parallel bilinear form for the finite element space
//     K = std::make_shared<mfem::ParBilinearForm>(fespace.get());

//     // Add a domain integrator for the diffusion term using the given diffusivity coefficient (cDx)
//     K->AddDomainIntegrator(new mfem::DiffusionIntegrator(*cDx));

//     // Assemble the bilinear form into a sparse matrix
//     K->Assemble();

//     // Temporary matrix to hold the assembled stiffness matrix
//     mfem::HypreParMatrix Khpm;

//     // Form the linear system, considering boundary conditions, parallel grid function (concentration, voxel), and the force term (F)
//     K->FormLinearSystem(boundary, parGF, F, Khpm, X1v, RHS);

//     // Convert the assembled stiffness matrix into a shared pointer for further use
//     Kmat = std::make_shared<mfem::HypreParMatrix>(Khpm);

//     // Scale the right-hand side vector for time-stepping by multiplying with the time step constant (Constants::dt)
//     RHS *= Constants::dt;

//     // Clean up dynamically allocated objects
//     // delete K;

// };


// void SolverSteps::ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, mfem::Array<int> &boundary, mfem::ProductCoefficient &m) 
// {
//     B = std::make_unique<mfem::ParLinearForm>(fespace.get());    

//     static mfem::ParGridFunction temp_parGF(fespace.get());
//     temp_parGF = parGF;

//     mfem::GridFunctionCoefficient coef(&temp_parGF);
//     B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

//     // Apply boundary conditions
//     B->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(m), boundary);

//     B->Assemble();
//     F = *B;
// }

// void SolverSteps::ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F) 
// {
//     B = std::make_unique<mfem::ParLinearForm>(fespace.get());    

//     static mfem::ParGridFunction temp_parGF(fespace.get());
//     temp_parGF = parGF;

//     mfem::GridFunctionCoefficient coef(&temp_parGF);
//     B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

//     B->Assemble();
//     F = *B;
// }

// void SolverSteps::SolverConditions(std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &solver, mfem::HypreSmoother &smoother){

//     // Set up the solver for the mass matrix.
//     solver.iterative_mode = false; // Use direct solving for the system matrix
//     solver.SetRelTol(1e-7); // Set relative tolerance for the solver
//     solver.SetAbsTol(0); // Set absolute tolerance for the solver
//     solver.SetMaxIter(102); // Limit the maximum number of iterations
//     solver.SetPrintLevel(0); // Suppress output from the solver
//     smoother.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
//     solver.SetPreconditioner(smoother); // Attach the preconditioner to the solver
//     solver.SetOperator(*Mmat); // Set the mass matrix as the operator to solve

// };

void SolverSteps::MassMatrix(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat) {
    mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
    if (!local_fespace) {
        throw std::runtime_error("fespace is null");
    }
    // Create a parallel bilinear form for the finite element space
    M = new mfem::ParBilinearForm(local_fespace);

    // Copy the provided potential field to a new grid function for internal operations
    temp_ps = new mfem::ParGridFunction(local_fespace);
    *temp_ps = psx;

    // Wrap the grid function in a coefficient (psi value) to use in the mass matrix integrator
    coef = new mfem::GridFunctionCoefficient(temp_ps);

    // Add a domain integrator to the bilinear form using the weighted mass integrator
    M->AddDomainIntegrator(new mfem::MassIntegrator(*coef));
    
    // Assemble the bilinear form into a sparse matrix representation
    M->Assemble();

    // Construct the mass matrix and store it in HPM
    mfem::HypreParMatrix HPM;
    M->FormSystemMatrix(boundary_dofs, HPM); // should be form linear system
    
    // Assign the mass matrix to the shared pointer
    Mmat = std::make_shared<mfem::HypreParMatrix>(HPM);

    // Clean up dynamically allocated objects
    delete coef; 
    delete temp_ps;
    delete M;
}

void SolverSteps::StiffnessMatrix(std::shared_ptr<mfem::GridFunctionCoefficient> cDx, mfem::Array<int> boundary, mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, std::shared_ptr<mfem::HypreParMatrix> &Kmat, mfem::HypreParVector &RHS){
    mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
    if (!local_fespace) {
        throw std::runtime_error("fespace is null");
    }
    // Create a parallel bilinear form for the finite element space
    K = std::make_shared<mfem::ParBilinearForm>(local_fespace);

    // Add a domain integrator for the diffusion term using the given diffusivity coefficient (cDx)
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(*cDx));

    // Assemble the bilinear form into a sparse matrix
    K->Assemble();

    // Temporary matrix to hold the assembled stiffness matrix
    mfem::HypreParMatrix Khpm;

    // Form the linear system, considering boundary conditions, parallel grid function (concentration, voxel), and the force term (F)
    K->FormLinearSystem(boundary, parGF, F, Khpm, X1v, RHS);

    // Convert the assembled stiffness matrix into a shared pointer for further use
    Kmat = std::make_shared<mfem::HypreParMatrix>(Khpm);

    // Scale the right-hand side vector for time-stepping by multiplying with the time step constant (Constants::dt)
    RHS *= Constants::dt;
}

void SolverSteps::ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, mfem::Array<int> &boundary, mfem::ProductCoefficient &m) {
    mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
    if (!local_fespace) {
        throw std::runtime_error("fespace is null");
    }
    B = std::make_unique<mfem::ParLinearForm>(local_fespace);    

    static mfem::ParGridFunction temp_parGF(local_fespace);
    temp_parGF = parGF;

    mfem::GridFunctionCoefficient coef(&temp_parGF);
    B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

    // Apply boundary conditions
    B->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(m), boundary);

    B->Assemble();
    F = *B;
}

void SolverSteps::ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F) {
    mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
    if (!local_fespace) {
        throw std::runtime_error("fespace is null");
    }
    B = std::make_unique<mfem::ParLinearForm>(local_fespace);    

    static mfem::ParGridFunction temp_parGF(local_fespace);
    temp_parGF = parGF;

    mfem::GridFunctionCoefficient coef(&temp_parGF);
    B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

    B->Assemble();
    F = *B;
}

void SolverSteps::SolverConditions(std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &solver, mfem::HypreSmoother &smoother){
    // Set up the solver for the mass matrix.
    solver.iterative_mode = false; // Use direct solving for the system matrix
    solver.SetRelTol(1e-8); // Set relative tolerance for the solver
    solver.SetAbsTol(0.0); // Set absolute tolerance for the solver
    solver.SetMaxIter(100); // Limit the maximum number of iterations
    solver.SetPrintLevel(0); // Suppress output from the solver
    smoother.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
    solver.SetPreconditioner(smoother); // Attach the preconditioner to the solver
    solver.SetOperator(*Mmat); // Set the mass matrix as the operator to solve
}