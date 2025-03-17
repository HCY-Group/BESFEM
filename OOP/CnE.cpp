/**
 * @file CnP.cpp
 * @brief Implementation of the electrolyte concentration class for battery simulations.
 */

#include "CnE.hpp"
#include "Constants.hpp"
#include "mfem.hpp"


CnE::CnE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), nbc_w_bdr(geo.nbc_w_bdr)
    
    {

    RxE = std::make_unique<mfem::ParGridFunction>(fespace.get());
    PeR = std::make_unique<mfem::ParGridFunction>(fespace.get());

    Mmate = std::make_shared<mfem::HypreParMatrix>();
    Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
    Me_prec.SetType(mfem::HypreSmoother::Jacobi);

    Kmate = std::make_shared<mfem::HypreParMatrix>();

    CeV0 = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    RHCe = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
    CeVn = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));

    }

void CnE::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
{
    Concentrations::SetInitialConcentration(Cn, initial_value);
    Concentrations::SetUpSolver(psx, Mmate, *Me_solver, Me_prec);

    // Apply Neumann boundary conditions to the reaction potential field
    ImposeNeumannBC(psx, *PeR);

}

void CnE::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{
    
    // Scale the reaction field by the transference number
    Concentrations::CreateReaction(Rx, *RxE, (-1.0 * Constants::t_minus));
    Concentrations::TotalReaction(*RxE, eCrnt);

    // Set up coefficients for Neumann boundary conditions
    mfem::ConstantCoefficient nbcCoef(infx); 
    mfem::GridFunctionCoefficient matCoef_R(PeR.get());
    mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);
    
    // Assemble the force term with boundary conditions applied
    Concentrations::ForceTerm(*RxE, ftE, nbc_w_bdr, m_nbcCoef, true); // true since applying BCs
    
    // Compute the diffusivity coefficient and stiffness matrix
    std::shared_ptr<mfem::GridFunctionCoefficient> cDe = Concentrations::Diffusivity(psx, Cn, false); // false using other equation
    eKx2 = std::make_shared<mfem::ParBilinearForm>(fespace.get());

    eKx2->Update();
    Concentrations::KMatrix(eKx2, boundary_dofs, Cn, ftE, Kmate, X1v, Feb, cDe);

    TmatR.reset(Add(1.0, *Mmate, -0.5 * Constants::dt, *Kmate));
    TmatL.reset(Add(1.0, *Mmate,  0.5 * Constants::dt, *Kmate));
 
    // Solve for the next time step concentration
    Cn.GetTrueDofs(*CeV0);	
    TmatR->Mult(*CeV0, *RHCe);
    *RHCe += Feb;
    
    Me_solver->SetOperator(*TmatL);
    Me_solver->Mult(*RHCe, *CeVn) ;

	Cn.Distribute(CeVn.get()); 
    
    // Update the total electrolyte salt conservation
    Concentrations::SaltConservation(Cn, psx);	
}
