

#include "../include/CnE.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include <optional>


CnE::CnE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), nbc_bdr(geo.nbc_bdr),
      De(fespace.get()), Rxe(fespace.get()), PeR(fespace.get()), cDe(&De), cAe(&Rxe), matCoef_R(&PeR),
      nbcCoef(0.0), Fet(fespace.get()), Me_solver(MPI_COMM_WORLD),
      CeV0(fespace.get()), RHCe(fespace.get()), CeVn(fespace.get()),
      Feb(fespace.get()), X1v(fespace.get())
    
    {

    }


void CnE::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    Concentrations::SetInitialConcentration(Cn, initial_value);

    // Assemble mass matrix
    mfem::GridFunctionCoefficient coef(&psx);
    SolverSteps::InitializeMassMatrix(coef, Me_init); 
    SolverSteps::FormSystemMatrix(Me_init, boundary_dofs, Mmate); 
    
    // Configure solver
    Me_solver.iterative_mode = false; // Enable iterative mode for the solver
    Me_prec.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
    SolverSteps::SolverConditions(Me_solver, Me_prec); // Set up the solver conditions for the mass matrix

    // Initialize stiffness operator (diffusivity)
    SolverSteps::InitializeStiffnessMatrix(cDe, Ke2); // Initialize

    // // Build boundary coefficient (Neumann BC weighting) HALF
    // PeR = psx;
    // PeR.Neg();
    // m_nbcCoef = std::make_unique<mfem::ProductCoefficient>(matCoef_R, nbcCoef);

    // Build force term (domain + boundary contributions)
    Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
    Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
    // Be_init->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*m_nbcCoef), nbc_bdr); // HALF
    Be_init->Assemble(); 
    Fet = *Be_init;

    boundary_dofs.SetSize(0);

    SolverSteps::FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the reaction potential
}

void CnE::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        // Assemble reaction source term
        Concentrations::CreateReaction(Rx, Rxe, (-1.0 * Constants::t_minus));
		cAe.SetGridFunction(&Rxe);
        Concentrations::TotalReaction(Rxe, eCrnt);

        nbcCoef.constant = infx;
		mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

        Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
        Be_init->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(m_nbcCoef), nbc_bdr);

        Be_init->Assemble();
        Fet = *Be_init; // Assign the updated linear form to Fet

        // salt diffusivity in the electrolyte					
		for (int vi = 0; vi < nV; vi++){
			// appendix equation A-21
			De(vi) = psx(vi)*Constants::D0*exp(-7.02-830*Cn(vi)+50000*Cn(vi)*Cn(vi));
		}
		cDe.SetGridFunction(&De);

        // Update stiffness operator
        SolverSteps::Update(Ke2); // Update the stiffness matrix with the new diffusivity coefficient
        SolverSteps::FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the updated values
   		Feb *= Constants::dt;

        // Crank-Nicolson matrices	
	    TmatR.reset(Add(1.0, Mmate, -0.5 * Constants::dt, Kmate));
        TmatL.reset(Add(1.0, Mmate,  0.5 * Constants::dt, Kmate));
 
        // Solve the system T·CeVn = RHCe
        Cn.GetTrueDofs(CeV0);
        TmatR->Mult(CeV0, RHCe);
        RHCe += Feb; // Add the right-hand side vector to the system

        Me_solver.SetOperator(*TmatL);
        Me_solver.Mult(RHCe, CeVn) ;

        // Recover updated concentration into GridFunction
	    Cn.Distribute(CeVn); 
        
    }	


void CnE::TimeStep(mfem::ParGridFunction &RxC, mfem::ParGridFunction &RxA, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        // Assemble reaction source term
        Concentrations::CreateReaction(RxC, RxA, Rxe, (-1.0 * Constants::t_minus));
		cAe.SetGridFunction(&Rxe);

        Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));

        Be_init->Assemble();
        Fet = *Be_init; // Assign the updated linear form to Fet

        // salt diffusivity in the electrolyte					
		for (int vi = 0; vi < nV; vi++){
			// appendix equation A-21
			De(vi) = psx(vi)*Constants::D0*exp(-7.02-830*Cn(vi)+50000*Cn(vi)*Cn(vi));
		}
		cDe.SetGridFunction(&De);

        // Update stiffness operator
        SolverSteps::Update(Ke2); // Update the stiffness matrix with the new diffusivity coefficient
        SolverSteps::FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the updated values
   		Feb *= Constants::dt;

        // Crank-Nicolson matrices	
	    TmatR.reset(Add(1.0, Mmate, -0.5 * Constants::dt, Kmate));
        TmatL.reset(Add(1.0, Mmate,  0.5 * Constants::dt, Kmate));
 
        // Solve the system T·CeVn = RHCe
        Cn.GetTrueDofs(CeV0);
        TmatR->Mult(CeV0, RHCe);
        RHCe += Feb; // Add the right-hand side vector to the system

        Me_solver.SetOperator(*TmatL);
        Me_solver.Mult(RHCe, CeVn) ;

        // Recover updated concentration into GridFunction
	    Cn.Distribute(CeVn); 
        
    }	