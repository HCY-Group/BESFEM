#include "../include/CnE.hpp"
#include "../include/Constants.hpp"
#include "../include/SimTypes.hpp"
#include "mfem.hpp"
#include <optional>


CnE::CnE(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::CellMode mode)
    : ConcentrationBase(geo, para), geometry(geo), domain_parameters(para), utils(geo, para), boundary_conditions(bc), fespace(geo.parfespace), nbc_bdr(bc.nbc_bdr), fem(geo.parfespace),
      De(fespace.get()), Rxe(fespace.get()), PeR(fespace.get()), cDe(&De), cAe(&Rxe), matCoef_R(&PeR),
      nbcCoef(0.0), Fet(fespace.get()), Me_solver(MPI_COMM_WORLD),
      CeV0(fespace.get()), RHCe(fespace.get()), CeVn(fespace.get()),
      Feb(fespace.get()), X1v(fespace.get()), mode_(mode), PsVc(fespace.get())
    
    {}


void CnE::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{

    if(mode_ == sim::CellMode::HALF){
       if (mfem::Mpi::WorldRank() == 0) {std::cout << "Initializing CnE for Half Cell" << std::endl;}

        utils.SetInitialValue(Cn, initial_value);

        // Assemble mass matrix
        mfem::GridFunctionCoefficient coef(&psx);
        fem.InitializeMassMatrix(coef, Me_init); 
        fem.FormSystemMatrix(Me_init, boundary_dofs, Mmate); 
        
        // Configure solver
        Me_solver.iterative_mode = false;
        Me_prec.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
        fem.SolverConditions(Me_solver, Me_prec); // Set up the solver conditions for the mass matrix

        // Initialize stiffness operator (diffusivity)
        fem.InitializeStiffnessMatrix(cDe, Ke2); // Initialize

        // Build boundary coefficient (Neumann BC weighting)
        PeR = psx;
        PeR.Neg();
        m_nbcCoef = std::make_unique<mfem::ProductCoefficient>(matCoef_R, nbcCoef);

        // Build force term (domain + boundary contributions)
        Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
        Be_init->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*m_nbcCoef), nbc_bdr);
        Be_init->Assemble();
        Fet = *Be_init;

        fem.FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the reaction potential
    
        // psx.GetTrueDofs(PsVc); // Extract true degrees of freedom
    }

    else if(mode_ == sim::CellMode::FULL){
        if (mfem::Mpi::WorldRank() == 0) {std::cout << "Initializing CnE for Full Cell" << std::endl;}

        utils.SetInitialValue(Cn, initial_value);

        // Assemble mass matrix
        mfem::GridFunctionCoefficient coef(&psx);
        fem.InitializeMassMatrix(coef, Me_init); 
        fem.FormSystemMatrix(Me_init, boundary_dofs, Mmate); 
        
        // Configure solver
        Me_solver.iterative_mode = false; // Enable iterative mode for the solver
        Me_prec.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
        fem.SolverConditions(Me_solver, Me_prec); // Set up the solver conditions for the mass matrix

        // Initialize stiffness operator (diffusivity)
        fem.InitializeStiffnessMatrix(cDe, Ke2); // Initialize

        // Build force term (domain + boundary contributions)
        Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
        Be_init->Assemble(); 
        Fet = *Be_init;

        fem.FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the reaction potential
    
        // psx.GetTrueDofs(PsVc); // Extract true degrees of freedom
    }
}

void CnE::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        // Assemble reaction source term
        utils.InitializeReaction(Rx, Rxe, (-1.0 * Constants::t_minus));
		cAe.SetGridFunction(&Rxe);
        utils.CalculateReactionInfx(Rxe, eCrnt);

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
        fem.Update(Ke2); // Update the stiffness matrix with the new diffusivity coefficient
        fem.FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the updated values
   		Feb *= Constants::dt;

        // Crank-Nicolson matrices	
	    TmatR.reset(Add(1.0, Mmate, -0.5 * Constants::dt, Kmate));
        TmatL.reset(Add(1.0, Mmate,  0.5 * Constants::dt, Kmate));
 
        // Solve the system T·CeVn = RHCe
        Cn.GetTrueDofs(CeV0);
        TmatR->Mult(CeV0, RHCe);
        RHCe += Feb; // Add the right-hand side vector to the system

        Me_solver.SetOperator(*TmatL);
        Me_solver.SetPreconditioner(Me_prec);
        Me_solver.Mult(RHCe, CeVn) ;

        // Clamp Negatives
        for (int p = 0; p < CeVn.Size(); p++){
            if (CeVn(p) < 0.0){
                (CeVn)(p) = 1.0e-6;} // prevent negative concentrations
        }

        

        // Recover updated concentration into GridFunction
	    Cn.Distribute(CeVn); 
        
    }	


void CnE::UpdateConcentration(mfem::ParGridFunction &RxC, mfem::ParGridFunction &RxA, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        // Assemble reaction source term
        utils.InitializeReaction(RxA, RxC, Rxe, (-1.0 * Constants::t_minus));
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
        // fem.Update(Ke2); // Update the stiffness matrix with the new diffusivity coefficient
        Ke2->Update();
        Ke2->Assemble();
        fem.FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); // Form the linear system for the updated values
   		Feb *= Constants::dt;

        // Crank-Nicolson matrices	
	    TmatR.reset(Add(1.0, Mmate, -0.5 * Constants::dt, Kmate));
        TmatL.reset(Add(1.0, Mmate,  0.5 * Constants::dt, Kmate));
 
        // Solve the system T·CeVn = RHCe
        Cn.GetTrueDofs(CeV0);
        TmatR->Mult(CeV0, RHCe);
        RHCe += Feb; // Add the right-hand side vector to the system

        Me_solver.SetOperator(*TmatL);
        // Me_solver.SetPreconditioner(Me_prec);

        Me_solver.Mult(RHCe, CeV0) ;

        // // Update only the solid region MAKE INTO FUNCTION
        // for (int p = 0; p < CeV0.Size(); p++){
        //     if (PsVc(p) < 3.0e-1){ // 1e-1 works, but gaps still get smaller
        //         (CeV0)(p) = Constants::init_CnE;} // Cp0 initial value
        // }

        // Recover updated concentration into GridFunction
	    Cn.Distribute(CeV0); 
    }	


void CnE::SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {

    CeC = 0.0; // Initialize total salt concentration to zero
    
    // Temporary grid function to hold the product of concentration and potential fields
   static mfem::ParGridFunction CeT(fespace.get());

    // Calculate the product of concentration and potential for each element
    CeT = Cn;
    CeT *= psx;
    
    // Loop through all elements to calculate the average concentration for each element
    for (int ei = 0; ei < nE; ei++){
        CeT.GetNodalValues(ei,VtxVal); // Retrieve the nodal values for the current element 
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}        
        EAvg(ei) = val/nC; // Compute the average concentration for the element
        // Update the total salt concentration by adding the weighted average for the element
        CeC += EAvg(ei)*EVol(ei); 
    }
    
    // Perform a global reduction to sum up the salt concentration across all MPI processes
    MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
    
    // Compute the average concentration across the entire electrolyte
    CeAvg = gCeC/gtPse;	

    // // Adjust the concentration field by normalizing the average concentration with the initial value

    // double damp = 1e-2;

    // for (int i = 0; i < Cn.Size(); i++) {
    //     Cn(i) -= damp * (CeAvg - Ce0) * psx(i); 
    // }    

    Cn -= (CeAvg-Ce0);
    MPI_Barrier(MPI_COMM_WORLD);

}