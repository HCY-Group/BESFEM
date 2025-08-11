#include "../include/SolverSteps.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

// Constructor for shared_ptr
SolverSteps::SolverSteps(std::shared_ptr<mfem::ParFiniteElementSpace> fespace)
: fespace(fespace), raw_fespace(nullptr), local_fespace(fespace.get()) {}

// Constructor for raw pointer
SolverSteps::SolverSteps(mfem::ParFiniteElementSpace* fespace)
: raw_fespace(fespace), fespace(nullptr), local_fespace(fespace) {}

// Destructor
SolverSteps::~SolverSteps() {}


void SolverSteps::InitializeMassMatrix(mfem::Coefficient &coef, std::unique_ptr<mfem::ParBilinearForm> &M) {

    M = std::make_unique<mfem::ParBilinearForm>(local_fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(coef));
    M->Assemble();

    std::cout << "Mass matrix initialized." << std::endl;

}

void SolverSteps::InitializeStiffnessMatrix(mfem::Coefficient &coef, std::unique_ptr<mfem::ParBilinearForm> &K) {

    K = std::make_unique<mfem::ParBilinearForm>(local_fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(coef));
    K->Assemble();

    std::cout << "Stiffness matrix initialized." << std::endl;

}

// void SolverSteps::InitializeForceTerm(mfem::Coefficient &coef, std::unique_ptr<mfem::ParLinearForm> &B) {

//     B = std::make_unique<mfem::ParLinearForm>(local_fespace);
//     B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));
//     B->Assemble();

//     std::cout << "Free energy initialized." << std::endl;
// }

void SolverSteps::InitializeForceTerm(
    mfem::Coefficient &coef,
    std::unique_ptr<mfem::ParLinearForm> &B,
    mfem::Coefficient *boundary_coef,
    mfem::Array<int> *boundary_attr)
{
    B = std::make_unique<mfem::ParLinearForm>(local_fespace);

    // Domain term
    B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

    // Optional boundary term
    if (boundary_coef && boundary_attr)
    {
        B->AddBoundaryIntegrator(
            new mfem::BoundaryLFIntegrator(*boundary_coef),
            *boundary_attr
        );
    }

    B->Assemble();

    // std::cout << "Force term initialized." << std::endl;
}


void SolverSteps::FormSystemMatrix(std::unique_ptr<mfem::ParBilinearForm> &M, mfem::Array<int> &boundary_dofs, std::shared_ptr<mfem::HypreParMatrix> &MassMatrix) {
    if (!M) {
        throw std::runtime_error("Bilinear form M is not initialized.");
    }
    
    // Temporary matrix to hold the assembled mass matrix
    mfem::HypreParMatrix A;

    // Form the system matrix from the bilinear form
    M->FormSystemMatrix(boundary_dofs, A);
    MassMatrix = std::make_shared<mfem::HypreParMatrix>(A);
    
    std::cout << "System matrix formed." << std::endl;
}

void SolverSteps::FormLinearSystem(std::unique_ptr<mfem::ParBilinearForm> &K, mfem::Array<int> &boundary_dofs, mfem::ParGridFunction &x, mfem::ParLinearForm &b, std::shared_ptr<mfem::HypreParMatrix> &StiffnessMatrix, mfem::HypreParVector &X, mfem::HypreParVector &B) {
    if (!K) {
        throw std::runtime_error("Bilinear form K is not initialized.");
    }

    // Temporary matrix to hold the assembled stiffness matrix
    mfem::HypreParMatrix A;

    K->FormLinearSystem(boundary_dofs, x, b, A, X, B);
    StiffnessMatrix = std::make_shared<mfem::HypreParMatrix>(A);

    // std::cout << "Linear system formed." << std::endl;
}

void SolverSteps::Update(std::unique_ptr<mfem::ParLinearForm> &B) {
    if (!B) {
        throw std::runtime_error("Linear form is not initialized.");
    }

    B->Update(); // Update the linear form with the current grid function values
    // Assemble the linear form
    B->Assemble();
}

void SolverSteps::Update(std::unique_ptr<mfem::ParBilinearForm> &B) {
    if (!B) {
        throw std::runtime_error("Bilinear form is not initialized.");
    }

    B->Update(); // Update the bilinear form with the current grid function values
    // Assemble the bilinear form
    B->Assemble();
}

void SolverSteps::SolverConditions(std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &solver, mfem::Solver &preconditioner){
    // Set up the solver for the mass matrix.
    solver.iterative_mode = false; // Use direct solving for the system matrix
    solver.SetRelTol(1e-7); // Set relative tolerance for the solver
    solver.SetAbsTol(0.0); // Set absolute tolerance for the solver
    solver.SetMaxIter(102); // Limit the maximum number of iterations
    solver.SetPrintLevel(0); // Suppress output from the solver
    solver.SetPreconditioner(preconditioner); // Attach the preconditioner to the solver
    solver.SetOperator(*Mmat); // Set the mass matrix as the operator to solve
}

void SolverSteps::SolverConditions(mfem::CGSolver &solver, mfem::Solver &preconditioner){
    // Set up the solver for the mass matrix.
    solver.iterative_mode = false; // Use direct solving for the system matrix
    solver.SetRelTol(1e-7); // Set relative tolerance for the solver
    solver.SetAbsTol(0.0); // Set absolute tolerance for the solver
    solver.SetMaxIter(102); // Limit the maximum number of iterations
    solver.SetPrintLevel(0); // Suppress output from the solver
    solver.SetPreconditioner(preconditioner); // Attach the preconditioner to the solver
}

// void SolverSteps::MassMatrix(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat) {
//     mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
//     if (!local_fespace) {
//         throw std::runtime_error("fespace is null");
//     }
//     // Create a parallel bilinear form for the finite element space
//     M = new mfem::ParBilinearForm(local_fespace);

//     // Copy the provided potential field to a new grid function for internal operations
//     temp_ps = new mfem::ParGridFunction(local_fespace);
//     *temp_ps = psx;

//     // Wrap the grid function in a coefficient (psi value) to use in the mass matrix integrator
//     coef = new mfem::GridFunctionCoefficient(temp_ps);

//     // Add a domain integrator to the bilinear form using the weighted mass integrator
//     M->AddDomainIntegrator(new mfem::MassIntegrator(*coef));
    
//     // Assemble the bilinear form into a sparse matrix representation
//     M->Assemble();

//     mfem::Array<int> boundary_dofs;
//     boundary_dofs.SetSize(0); 

//     // Construct the mass matrix and store it in HPM
//     mfem::HypreParMatrix HPM;
//     M->FormSystemMatrix(boundary_dofs, HPM); // should be form linear system
    
//     // Assign the mass matrix to the shared pointer
//     Mmat = std::make_shared<mfem::HypreParMatrix>(HPM);

//     // Clean up dynamically allocated objects
//     delete coef; 
//     delete temp_ps;
//     delete M;
// }


// void SolverSteps::StiffnessMatrix(std::shared_ptr<mfem::GridFunctionCoefficient> cDx, mfem::Array<int> boundary, mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, std::shared_ptr<mfem::HypreParMatrix> &Kmat, mfem::HypreParVector &RHS){
//     mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
//     if (!local_fespace) {
//         throw std::runtime_error("fespace is null");
//     }
//     // Create a parallel bilinear form for the finite element space
//     K = std::make_shared<mfem::ParBilinearForm>(local_fespace);

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
// }


// void SolverSteps::StiffnessMatrix(mfem::Coefficient &cDx, mfem::Array<int> boundary, std::shared_ptr<mfem::HypreParMatrix> &Kmat){
//     mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
//     if (!local_fespace) {
//         throw std::runtime_error("fespace is null");
//     }
//     // Create a parallel bilinear form for the finite element space
//     K = std::make_shared<mfem::ParBilinearForm>(local_fespace);

//     // Add a domain integrator for the diffusion term using the given diffusivity coefficient (cDx)
//     K->AddDomainIntegrator(new mfem::DiffusionIntegrator(cDx));

//     // Assemble the bilinear form into a sparse matrix
//     K->Assemble();

//     // Temporary matrix to hold the assembled stiffness matrix
//     mfem::HypreParMatrix Khpm;

//     // Dummy solution and RHS for Form Linear System
//     mfem::GridFunction dummy_sol(local_fespace);
//     mfem::Vector dummy_rhs, dummy_X;

//     dummy_sol = 0.0;

//     K->FormLinearSystem(boundary, dummy_sol, dummy_sol, Khpm, dummy_X, dummy_rhs);

//     // Convert the assembled stiffness matrix into a shared pointer for further use
//     Kmat = std::make_shared<mfem::HypreParMatrix>(Khpm);
// }

// void SolverSteps::ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, mfem::Array<int> &boundary, mfem::ProductCoefficient &m) {
//     mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
//     if (!local_fespace) {
//         throw std::runtime_error("fespace is null");
//     }
//     B = std::make_unique<mfem::ParLinearForm>(local_fespace);    

//     static mfem::ParGridFunction temp_parGF(local_fespace);
//     temp_parGF = parGF;

//     mfem::GridFunctionCoefficient coef(&temp_parGF);
//     B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

//     // Apply boundary conditions
//     B->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(m), boundary);

//     B->Assemble();
//     F = *B;
// }

// void SolverSteps::ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F) {
//     mfem::ParFiniteElementSpace* local_fespace = fespace ? fespace.get() : raw_fespace;
//     if (!local_fespace) {
//         throw std::runtime_error("fespace is null");
//     }
//     B = std::make_unique<mfem::ParLinearForm>(local_fespace);    

//     static mfem::ParGridFunction temp_parGF(local_fespace);
//     temp_parGF = parGF;

//     mfem::GridFunctionCoefficient coef(&temp_parGF);
//     B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

//     B->Assemble();
//     F = *B;
// }

