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

    // std::cout << "Mass matrix initialized." << std::endl;

}

void SolverSteps::InitializeStiffnessMatrix(mfem::Coefficient &coef, std::unique_ptr<mfem::ParBilinearForm> &K) {

    K = std::make_unique<mfem::ParBilinearForm>(local_fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(coef));
    K->Assemble();

    // std::cout << "Stiffness matrix initialized." << std::endl;

}

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


void SolverSteps::FormSystemMatrix(std::unique_ptr<mfem::ParBilinearForm> &M, mfem::Array<int> &boundary_dofs, mfem::HypreParMatrix &A) {
    if (!M) {
        throw std::runtime_error("Bilinear form M is not initialized.");
    }
    
    M->FormSystemMatrix(boundary_dofs, A);

}

void SolverSteps::FormLinearSystem(std::unique_ptr<mfem::ParBilinearForm> &K, mfem::Array<int> &boundary_dofs, mfem::ParGridFunction &x, mfem::ParLinearForm &b, mfem::HypreParMatrix &A, mfem::HypreParVector &X, mfem::HypreParVector &B) {
    if (!K) {
        throw std::runtime_error("Bilinear form K is not initialized.");
    }

    K->FormLinearSystem(boundary_dofs, x, b, A, X, B);
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

void SolverSteps::SolverConditions(mfem::HypreParMatrix &Mmat, mfem::CGSolver &solver, mfem::Solver &preconditioner){
    // Set up the solver for the mass matrix.
    solver.iterative_mode = false; // Use direct solving for the system matrix
    solver.SetRelTol(1e-6); // Set relative tolerance for the solver
    solver.SetAbsTol(0.0); // Set absolute tolerance for the solver
    solver.SetMaxIter(102); // Limit the maximum number of iterations
    solver.SetPrintLevel(0); // Suppress output from the solver
    solver.SetPreconditioner(preconditioner); // Attach the preconditioner to the solver
    solver.SetOperator(Mmat); // Set the mass matrix as the operator to solve
}

void SolverSteps::SolverConditions(mfem::CGSolver &solver, mfem::Solver &preconditioner){
    // Set up the solver for the mass matrix.
    solver.iterative_mode = false; // Use direct solving for the system matrix
    solver.SetRelTol(1e-6); // Set relative tolerance for the solver
    solver.SetAbsTol(0.0); // Set absolute tolerance for the solver
    solver.SetMaxIter(102); // Limit the maximum number of iterations
    solver.SetPrintLevel(0); // Suppress output from the solver
    solver.SetPreconditioner(preconditioner); // Attach the preconditioner to the solver
}
