/**
 * @file CnP.cpp
 * @brief Implementation of the electrolyte concentration class for battery simulations.
 */

#include "CnE.hpp"
#include "../code/Constants.hpp"
#include "mfem.hpp"
#include <optional>



CnE::CnE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), nbc_w_bdr(geo.nbc_w_bdr), Rxe(fespace.get()), De(fespace.get())
    
    {

    PeR = std::make_unique<mfem::ParGridFunction>(fespace.get());

    Me_init = std::make_unique<mfem::ParBilinearForm>(fespace.get());
    Mmate = std::make_shared<mfem::HypreParMatrix>();
    Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Me_prec.SetType(mfem::HypreSmoother::Jacobi);

    Kmate = std::make_shared<mfem::HypreParMatrix>();

    CeV0 = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    RHCe = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    CeVn = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));

    Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());

    De = mfem::ParGridFunction(fespace.get());
    Rxe = mfem::ParGridFunction(fespace.get());
    Fet = mfem::ParLinearForm(fespace.get());
    Feb = mfem::HypreParVector(fespace.get());

    cAe = mfem::GridFunctionCoefficient(&Rxe);
    cDe = mfem::GridFunctionCoefficient(&De);

    }

void CnE::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    Concentrations::SetInitialConcentration(Cn, initial_value);

    mfem::GridFunctionCoefficient coef(&psx);
    SolverSteps::InitializeMassMatrix(coef, Me_init); 
    SolverSteps::FormSystemMatrix(Me_init, boundary_dofs, Mmate); 
    SolverSteps::SolverConditions(Mmate, *Me_solver, Me_prec); // Set up the solver conditions for the mass matrix

    SolverSteps::InitializeStiffnessMatrix(cDe, Ke2); // Initialize

    Concentrations::ImposeNeumannBC(psx, *PeR); // Apply Neumann boundary conditions to the reaction potential field

    // Set up coefficients for Neumann boundary conditions
    mfem::ConstantCoefficient nbcCoef(infx); 
    mfem::GridFunctionCoefficient matCoef_R(PeR.get());
    mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    SolverSteps::InitializeForceTerm(cAe, Be_init); // Initialize the force term for the reaction potential
    Fet = std::move(*Be_init); // Move the initialized force term to Fet

    SolverSteps::FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the reaction potential

}

void CnE::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        Concentrations::CreateReaction(Rx, Rxe, (-1.0 * Constants::t_minus));
        cAe.SetGridFunction(&Rxe); // Set the reaction term coefficient for the force term

        Concentrations::TotalReaction(Rxe, eCrnt);

        // Set up coefficients for Neumann boundary conditions
        mfem::ConstantCoefficient nbcCoef(infx); 
        mfem::GridFunctionCoefficient matCoef_R(PeR.get());
        mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

        SolverSteps::Update(Be_init); // Update the force term with the current reaction term
        Fet = std::move(*Be_init);

        // Compute the diffusivity coefficient and stiffness matrix
        for (int vi = 0; vi < nV; vi++) {
            De(vi) = psx(vi) * Constants::D0 * exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
        }
        cDe.SetGridFunction(&De); // Set the diffusivity coefficient for the stiffness matrix with updated concentration
     
        SolverSteps::Update(Ke2); // Update the stiffness matrix with the new diffusivity coefficient
        SolverSteps::FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the updated values
        Feb *= Constants::dt; // Scale the right-hand side vector by the time step

        // Crank-Nicolson time-stepping
        TmatR.reset(Add(1.0, *Mmate, -0.5 * Constants::dt, *Kmate));
        TmatL.reset(Add(1.0, *Mmate,  0.5 * Constants::dt, *Kmate));
        
        // Solve the system T·CeVn = RHCe
        Cn.GetTrueDofs(*CeV0);
        TmatR->Mult(*CeV0, *RHCe);
        *RHCe += Feb; // Add the right-hand side vector to the system
        
        Me_solver->SetOperator(*TmatL);
        Me_solver->Mult(*RHCe, *CeVn) ;

	    Cn.Distribute(CeVn.get()); 
        
        // Update the total electrolyte salt conservation
        Concentrations::SaltConservation(Cn, psx);
    }	
        
        
