#include "Concentrations_Base.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"

/**
 * @class Concentrations
 * @brief Base class for managing concentration calculations in a finite element framework.
 * 
 * This class provides common functionalities for initializing, updating, and solving 
 * concentration-related equations in parallel simulations.
 */
Concentrations::Concentrations(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), EVol(mh.GetElementVolume()), gtPsi(mesh_handler.GetTotalPsi()), 
      gtPse(mesh_handler.GetTotalPse())

{
    
    nE = mesh_handler.GetNE(); ///< Number of elements in the mesh.
    nC = mesh_handler.GetNC(); ///< Number of nodes per element.
    nV = mesh_handler.GetNV(); ///< Total number of vertices in the mesh.

    // mfem::ParGridFunction CeT(fespace); // I think this is the issue

    // CeT = new mfem::ParGridFunction(fespace);


}

/**
 * @brief Sets the initial values for concentration and related parameters.
 * 
 * @param Cn Grid function for concentration to initialize.
 * @param initial_value Initial value for concentration.
 * @param psx Grid function for psi values.
 * @param perform_lithiation Flag to perform lithiation calculations.
 */
void Concentrations::SetInitialValues(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation) {
    
    // Shared initialization tasks
    SetInitialConcentration(Cn, initial_value);

    if (perform_lithiation) {
        LithiationCalculation(Cn, psx);
    }

}

/**
 * @brief Sets initial concentration values for a given grid function.
 * 
 * @param Cn Grid function for concentration.
 * @param initial_value Initial value to set.
 */
void Concentrations::SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value) {
    
    for (int i = 0; i < Cn.Size(); ++i) {
        Cn(i) = initial_value;
    }

}

/**
 * @brief Sets up the solver for concentration calculations.
 * 
 * @param psx Grid function for psi values.
 * @param Mmat Matrix to store the mass matrix.
 * @param m_solver Conjugate gradient solver.
 * @param smoother Preconditioner for the solver.
 */
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

/**
 * @brief Calculates lithiation-related updates for the given grid functions.
 * 
 * @param Cn Concentration grid function.
 * @param psx Psi grid function affecting lithiation.
 */
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

/**
 * @brief Applies Neumann boundary conditions to the grid function.
 * 
 * @param psx Input grid function for psi values.
 * @param PGF Output grid function for boundary conditions.
 */
void Concentrations::ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF) {
    PGF = psx;
    PGF.Neg();
}

/**
 * @brief Creates a reaction term by scaling a grid function.
 * 
 * @param Rx1 Input reaction rate grid function.
 * @param Rx2 Output reaction rate grid function.
 * @param value Scaling factor for the reaction.
 */
void Concentrations::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {

    Rx2 = Rx1;
    Rx2 *= value;

    // std::cout << "Rxn in Create Reaction in Base Class" << std::endl;
    // Rx2.Print(std::cout);

}

/**
 * @brief Computes a force term for the concentration system.
 * 
 * @param gfc Grid function for concentration.
 * @param Fxx Output linear form for force terms.
 * @param boundary Array of boundary conditions.
 * @param m Product coefficient for the force term.
 * @param apply_boundary_conditions Flag to apply boundary conditions.
 */
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

    // Bx2->Print(std::cout);

}

/**
 * @brief Computes the total reaction current.
 * 
 * This method integrates the reaction rate over all elements to calculate
 * the total reaction current. The result is adjusted based on the mesh's
 * physical dimensions.
 * 
 * @param Rx Grid function representing the reaction rate.
 * @param xCrnt Reference to the variable holding the computed reaction current.
 */
void Concentrations::TotalReaction(mfem::ParGridFunction &Rx, double xCrnt) {
    
    xCrnt = 0.0;
    Array<double> VtxVal(nC);
    Vector EAvg(nE);
    for (int ei = 0; ei < nE; ei++) {
        Rx.GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / nC;
        xCrnt += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&xCrnt, &geCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    infx = geCrnt / (mesh_handler.L_w);

    // std::cout << "geCrnt: " << geCrnt << std::endl;
    // std::cout << "infx: " << infx << std::endl;



}

/**
 * @brief Computes the diffusivity based on input parameters.
 * 
 * The diffusivity is computed using either a particle or electrolyte equation,
 * depending on the input flag. The results are stored in a grid function coefficient.
 * 
 * @param psx Grid function representing psi values.
 * @param Cn Grid function for concentration values.
 * @param particle_electrolyte Flag indicating whether to use the particle or electrolyte equation.
 * @return Shared pointer to the computed grid function coefficient for diffusivity.
 */
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

/**
 * @brief Constructs the stiffness matrix for the concentration system.
 * 
 * This method assembles the stiffness matrix using diffusion integrators
 * and forms the linear system for the solver.
 * 
 * @param boundary Array specifying boundary conditions.
 * @param Cn Grid function for concentration values.
 * @param Fxx Linear form for the force term.
 * @param Kmatx Shared pointer to the resulting stiffness matrix.
 * @param X1v Input vector for the linear system.
 * @param Fxb Right-hand side vector for the system.
 * @param cDx Coefficient representing diffusivity.
 */
void Concentrations::KMatrix(mfem::Array<int> boundary, mfem::ParGridFunction &Cn, mfem::ParLinearForm &Fxx, std::shared_ptr<mfem::HypreParMatrix> &Kmatx, mfem::HypreParVector &X1v, mfem::HypreParVector &Fxb, mfem::GridFunctionCoefficient *cDx) {

    std::unique_ptr<mfem::ParBilinearForm> Kx2(new mfem::ParBilinearForm(fespace));

    mfem::HypreParMatrix Khpm;
    
    Kx2->AddDomainIntegrator(new DiffusionIntegrator(*cDx));
    Kx2->Assemble();
    Kx2->FormLinearSystem(boundary, Cn, Fxx, Khpm, X1v, Fxb);

    Kmatx = std::make_shared<mfem::HypreParMatrix>(Khpm);

    Fxb *= Constants::dt;


}

/**
 * @brief Performs salt conservation calculations.
 * 
 * This method adjusts the concentration values to enforce conservation of salt
 * within the electrolyte by computing the average and modifying the values accordingly.
 * 
 * @param Cn Grid function representing the current concentration values.
 * @param psx Grid function for psi values.
 */
void Concentrations::SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {

    CeC = 0.0;
    mfem::ParGridFunction CeT(fespace);
    Array<double> VtxVal(nC);
    Vector EAvg(nE);


    CeT = Cn;
    CeT *= psx;

    for (int ei = 0; ei < nE; ei++){
        CeT.GetNodalValues(ei,VtxVal) ;
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
        EAvg(ei) = val/nC;	
        CeC += EAvg(ei)*EVol(ei) ;
    }

    MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
    
    // average CnE throughout electrolyte
    CeAvg = gCeC/gtPse;	
    
    // adjust CnE
    Cn -= (CeAvg-Ce0);
    MPI_Barrier(MPI_COMM_WORLD);



}

// void Concentrations::Save(mfem::ParGridFunction &gf, const std::string &base_name) {
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the MPI rank
// 
//     std::string file_name = base_name + "." + std::to_string(rank) + ".gf";  // Use the dynamic base name
// 
//     std::ofstream ofs(file_name.c_str());
//     if (ofs.is_open()) {
//         gf.Save(ofs);  
//         ofs.close();
//     } else {
//         mfem::mfem_error("Error opening file to save.");
//     }
// }
