#include "CnP.hpp"
#include "mfem.hpp"

/**
 * @class CnP
 * @brief Manages particle concentration calculations.
 *
 * This class is derived from `Concentrations` and includes functionality for:
 * - Setting initial values.
 * - Time-stepping methods.
 */
CnP::CnP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh)
    
    {

    PsVc = mfem::HypreParVector(fespace); 
    RxP = new ParGridFunction(fespace); ///< Grid function for reaction rate.

    Mmatp = std::make_shared<mfem::HypreParMatrix>(); ///< Mass matrix for particle.
    Mp_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD); ///< Conjugate gradient solver.
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi); ///< Preconditioner for CG solver.

    Kmatp = std::make_shared<mfem::HypreParMatrix>(); ///< Stiffness matrix for electrolyte.
    Fcb = HypreParVector(fespace);

    CpV0 = new mfem::HypreParVector(fespace); ///< Initial particle concentration vector.
    RHCp = new mfem::HypreParVector(fespace); ///< Right-hand side vector for concentration updates.
    CpVn = new mfem::HypreParVector(fespace); ///< Updated particle concentration vector.


    }
    
/**
 * @fn void CnP::Initialize(mfem::ParGridFunction&, double, mfem::ParGridFunction&, bool)
 * @brief Sets initial values for concentration and related parameters.
 * 
 * @param Cn The concentration grid function to initialize.
 * @param initial_value Initial concentration value.
 * @param psx The grid function for psi related to the particle concentration.
 * @param perform_lithiation If true, perform lithiation calculations.
 */
void CnP::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation)
{
    Concentrations::SetInitialValues(Cn, initial_value, psx, perform_lithiation);
    Concentrations::SetUpSolver(psx, Mmatp, *Mp_solver, Mp_prec);

    psx.GetTrueDofs(PsVc);
}

/**
 * @brief Performs a time step for particle concentration calculations.
 * 
 * Updates the concentration using Backward Euler, including reaction
 * and diffusion terms.
 * 
 * @param Rx Reaction rate grid function.
 * @param Cn Concentration grid function.
 * @param psx Psi grid function.
 */
void CnP::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
{

    Concentrations::CreateReaction(Rx, *RxP, (1/Constants::rho));

    // Create dummy values to use Force Term Function 
    Array<int> dummy_boundary;
    mfem::ConstantCoefficient coef1(0.0);
    mfem::ConstantCoefficient coef2(0.0);
    mfem::ProductCoefficient dummy_coef(coef1, coef2);

    Concentrations::ForceTerm(*RxP, ftPC, dummy_boundary, dummy_coef, false); // false since not applying BCs

    std::shared_ptr<GridFunctionCoefficient> cDp = Concentrations::Diffusivity(psx, Cn, true); // true since using first equation
    Concentrations::KMatrix(boundary_dofs, Cn, ftPC, Kmatp, X1v, Fcb, cDp.get());

    Tmatp = Add(1.0, *Mmatp, -(Constants::dt), *Kmatp);

    int nDof = CpV0->Size();
    Cn.GetTrueDofs(*CpV0);

    Tmatp->Mult(*CpV0, *RHCp);
    *RHCp += Fcb;

    Mp_solver->Mult(*RHCp, *CpVn);

    // Update only the solid region MAKE INTO FUNCTION
    for (int p = 0; p < nDof; p++){
        if (PsVc(p) < 1.0e-5){
            (*CpVn)(p) = 0.3;} // Cp0 initial value
    }

    Cn.Distribute(CpVn);


    // Degree of Lithiation
    Concentrations::LithiationCalculation(Cn, psx);


}

