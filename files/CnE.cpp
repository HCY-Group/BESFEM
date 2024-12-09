#include "CnE.hpp"
#include "mfem.hpp"

/**
 * @class CnE
 * @brief Manages electrolyte concentration calculations.
 * 
 * This class is derived from `Concentrations` and includes functionality for:
 * - Setting initial values.
 * - Time-stepping methods using Crank-Nicolson.
 */
CnE::CnE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh) 
    
    {

    PeR = new mfem::ParGridFunction(fespace); ///< Grid function for reaction potential.
    RxE = new mfem::ParGridFunction(fespace); ///< Grid function for reaction rate.

    Mmate = std::make_shared<mfem::HypreParMatrix>(); ///< Mass matrix for electrolyte.
    Me_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD); ///< Conjugate gradient solver.
    Me_prec.SetType(mfem::HypreSmoother::Jacobi); ///< Preconditioner for CG solver. 

    Kmate = std::make_shared<mfem::HypreParMatrix>(); ///< Stiffness matrix for electrolyte.

    CeV0 = new mfem::HypreParVector(fespace); ///< Initial electrolyte concentration vector.
    RHCe = new mfem::HypreParVector(fespace); ///< Right-hand side vector for concentration updates.
    CeVn = new mfem::HypreParVector(fespace); ///< Updated electrolyte concentration vector.
 
    }

mfem::HypreParVector* CnE::CeVn = nullptr; // static variable to be used in reaction

/**
 * @brief Initializes the electrolyte concentration and solver setup.
 * 
 * @param Cn Grid function representing the initial concentration.
 * @param initial_value The initial value to set for concentration.
 * @param psx The grid function for psi related to the electrolyte concentration.
 * @param perform_lithiation Flag indicating whether to perform lithiation calculations.
 */
void CnE::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation)
{
    Concentrations::SetInitialValues(Cn, initial_value, psx, perform_lithiation);
    Concentrations::SetUpSolver(psx, Mmate, *Me_solver, Me_prec);

    ImposeNeumannBC(psx, *PeR);

}

/**
 * @brief Performs a time step for electrolyte concentration calculations.
 * 
 * Updates the concentration using Crank-Nicolson integration, including reaction
 * and diffusion terms.
 * 
 * @param Rx Reaction rate grid function.
 * @param Cn Concentration grid function.
 * @param psx Psi grid function.
 */
void CnE::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{

    Concentrations::CreateReaction(Rx, *RxE, (-1.0 * Constants::t_minus));
    Concentrations::TotalReaction(*RxE, eCrnt);

    // Values used in Force Term Function
    mfem::ConstantCoefficient nbcCoef(infx); 
    mfem::GridFunctionCoefficient matCoef_R(PeR);
    mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    mfem::Array<int> nbc_w_bdr(pmesh->bdr_attributes.Max());
	nbc_w_bdr = 0; 
	nbc_w_bdr[0] = 1;

    Concentrations::ForceTerm(*RxE, ftE, nbc_w_bdr, m_nbcCoef, true); // true since applying BCs
    std::shared_ptr<GridFunctionCoefficient> cDe = Concentrations::Diffusivity(psx, Cn, false); // false using other equation
    Concentrations::KMatrix(boundary_dofs, Cn, ftE, Kmate, X1v, Feb, cDe.get());

    // Crank-Nicolson matrices
    TmatR = Add(1.0,*Mmate, -0.5*Constants::dt, *Kmate);		
    TmatL = Add(1.0, *Mmate,  0.5*Constants::dt, *Kmate);	
    
    Cn.GetTrueDofs(*CeV0);	

    TmatR->Mult(*CeV0, *RHCe);
    *RHCe += Feb;
    
    Me_solver->SetOperator(*TmatL);
    Me_solver->Mult(*RHCe, *CeVn) ;

	Cn.Distribute(CeVn); 

    Concentrations::SaltConservation(Cn, psx);	


}
