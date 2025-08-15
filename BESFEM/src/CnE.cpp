/**
 * @file CnP.cpp
 * @brief Implementation of the electrolyte concentration class for battery simulations.
 */

#include "../include/CnE.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <optional>



CnE::CnE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), nbc_w_bdr(geo.nbc_w_bdr), Rxe(fespace.get()), De(fespace.get()),
        nbcCoef(0.0), matCoef_R(nullptr), m_nbcCoef(nullptr)
    
    {

    PeR = std::make_unique<mfem::ParGridFunction>(fespace.get());

    Me_init = std::make_unique<mfem::ParBilinearForm>(fespace.get());
    Mmate = std::make_shared<mfem::HypreParMatrix>();
    Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Me_prec.SetType(mfem::HypreSmoother::Jacobi);

    Kmate = std::make_shared<mfem::HypreParMatrix>();
    Ke2 = std::make_unique<mfem::ParBilinearForm>(fespace.get());

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
    
    Me_prec.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
    SolverSteps::SolverConditions(*Me_solver, Me_prec); // Set up the solver conditions for the mass matrix
    // SolverSteps::SolverConditions(Mmate, *Me_solver, Me_prec); // Set up the solver conditions for the mass matrix

    SolverSteps::InitializeStiffnessMatrix(cDe, Ke2); // Initialize

    *PeR = psx;
    PeR->Neg();
    
    nbcCoef = mfem::ConstantCoefficient(0.0);  // initialize with zero
    matCoef_R = mfem::GridFunctionCoefficient(PeR.get());
    m_nbcCoef = std::make_unique<mfem::ProductCoefficient>(matCoef_R, nbcCoef);

    // Construct Be_init once with the correct integrators
    Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
    Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
    Be_init->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*m_nbcCoef), nbc_w_bdr);
    Be_init->Assemble();

    // SolverSteps::InitializeForceTerm(cAe, Be_init); // Initialize the force term for the reaction potential
    // SolverSteps::InitializeForceTerm(cAe, Be_init, m_nbcCoef.get(), &nbc_w_bdr); // Initialize the force term for the reaction potential
    Fet = *Be_init;

    SolverSteps::FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the reaction potential

}

void CnE::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        Concentrations::CreateReaction(Rx, Rxe, (-1.0 * Constants::t_minus));
		cAe.SetGridFunction(&Rxe);
        
        Concentrations::TotalReaction(Rxe, eCrnt);

        // !!!!!!!
        // remake boundary condition for Neumann CnE
        // is the Neumann BC being updated and used correctly in each timestep?

        nbcCoef = mfem::ConstantCoefficient(infx);
		mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

        SolverSteps::Update(Be_init); // Update the linear form with the new coefficients
        Fet = std::move(*Be_init); // Move the contents of Be_init into Fet

        // salt diffusivity in the electrolyte					
		for (int vi = 0; vi < nV; vi++){
			// appendix equation A-21
			De(vi) = psx(vi)*Constants::D0*exp(-7.02-830*Cn(vi)+50000*Cn(vi)*Cn(vi));
		}
		cDe.SetGridFunction(&De);

        SolverSteps::Update(Ke2); // Update the stiffness matrix with the new diffusivity coefficient
        SolverSteps::FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the updated values
   		Feb *= Constants::dt;

        // Crank-Nicolson matrices	
	    TmatR.reset(Add(1.0, *Mmate, -0.5 * Constants::dt, *Kmate));
        TmatL.reset(Add(1.0, *Mmate,  0.5 * Constants::dt, *Kmate));
 
        // Solve the system T·CeVn = RHCe
        Cn.GetTrueDofs(*CeV0);
        TmatR->Mult(*CeV0, *RHCe);
        *RHCe += Feb; // Add the right-hand side vector to the system

        Me_solver->SetOperator(*TmatL);
        Me_solver->Mult(*RHCe, *CeVn) ;

	    Cn.Distribute(CeVn.get()); 

        // // After updating CnE:
        // for (int vi = 0; vi < nV; vi++) {
        //     if (Cn(vi) < 0.0) Cn(vi) = 0.0;
        //     if (Cn(vi) > 1.0) Cn(vi) = 1.0;
        // }
        
    }	

