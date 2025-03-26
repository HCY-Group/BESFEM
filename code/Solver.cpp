#include "Solver.hpp"
#include "Initialize_Geometry.hpp"
#include "Constants.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

// Constructor
Solver::Solver(Initialize_Geometry &geo)
: pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo)

{}

// Destructor
Solver::~Solver() {}

void Solver::MassMatrix(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat) {

    // Create a parallel bilinear form for the finite element space
    M = new mfem::ParBilinearForm(fespace.get());

    // Copy the provided potential field to a new grid function for internal operations
    temp_ps = new mfem::ParGridFunction(fespace.get());
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

};

void Solver::StiffnessMatrix(std::shared_ptr<mfem::GridFunctionCoefficient> cDx, mfem::Array<int> boundary, mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, std::shared_ptr<mfem::HypreParMatrix> &Kmat, mfem::HypreParVector &RHS){

    // Create a parallel bilinear form for the finite element space
    K = std::make_shared<mfem::ParBilinearForm>(fespace.get());

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

    // Clean up dynamically allocated objects
    // delete K;

};

void Solver::ForceTerm(mfem::ParGridFunction &parGF, mfem::ParLinearForm &F, mfem::Array<int> boundary, mfem::ProductCoefficient m, bool apply_boundary_conditions){

    B = std::make_unique<mfem::ParLinearForm>(fespace.get());    

    // Initialize a new ParGridFunction to hold the input field and copy values from parGF
    static mfem::ParGridFunction temp_parGF(fespace.get());
    temp_parGF = parGF;
    
    // Create a GridFunctionCoefficient from the ParGridFunction for use in integrators
    mfem::GridFunctionCoefficient coef(&temp_parGF);

    // Add a domain integrator to compute contributions from the entire domain
    B->AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));

    // If boundary conditions are to be applied, add a boundary integrator
    if (apply_boundary_conditions) {
        B->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(m), boundary);
    }
    
    // Assemble the linear form to finalize its representation
    B->Assemble();
    
    // Move the assembled linear form into the provided force reference
    F = *B;

};

void Solver::SolverConditions(std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &solver, mfem::HypreSmoother &smoother){

    // Set up the solver for the mass matrix.
    solver.iterative_mode = false; // Use direct solving for the system matrix
    solver.SetRelTol(1e-7); // Set relative tolerance for the solver
    solver.SetAbsTol(0); // Set absolute tolerance for the solver
    solver.SetMaxIter(102); // Limit the maximum number of iterations
    solver.SetPrintLevel(0); // Suppress output from the solver
    smoother.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
    solver.SetPreconditioner(smoother); // Attach the preconditioner to the solver
    solver.SetOperator(*Mmat); // Set the mass matrix as the operator to solve

};

