#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"
#include "SimTypes.hpp"

// /**
//  * @file CnE.hpp
//  * @brief Electrolyte concentration (CnE) solver for battery simulations.
//  *
//  * Provides initialization and Crank–Nicolson time-stepping for the electrolyte
//  * concentration field using MFEM/Hypre operators, reaction coupling, and
//  * Neumann boundary conditions on a marked boundary set.
//  */

class Initialize_Geometry;
class Domain_Parameters;
class BoundaryConditions;
class CellMode;

// /**
//  * @defgroup concentrations Concentration Modules
//  * @brief Classes that advance concentration fields.
//  * @{
//  */


// /**
//  * @class CnE
//  * @brief Electrolyte concentration solver (electrolyte phase).
//  * @ingroup concentrations
//  *
//  * @details
//  * - Assembles mass/stiffness operators in parallel (MFEM/Hypre).
//  * - Handles domain and Neumann boundary forcing terms.
//  * - Uses a Crank–Nicolson update for time integration.
//  */
class CnE : public ConcentrationBase {
public:

    // /**
    //  * @brief Construct the electrolyte solver.
    //  * @param geo  Geometry/space container (mesh, FESpace, boundary markers).
    //  * @param para Domain/physics parameters (time step, constants, etc.).
    //  */
    CnE(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::CellMode mode);

    // /**
    //  * @brief Initialize solver and assemble operators.
    //  *
    //  * Sets the initial concentration, assembles mass and stiffness operators,
    //  * prepares boundary/domain source terms, and configures solver and preconditioner.
    //  *
    //  * @param Cn            Electrolyte concentration grid function (in/out).
    //  * @param initial_value Initial scalar value used for \p Cn.
    //  * @param psx           Phase field ψ used for masking/weights/BCs.
    //  */
    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    // /**
    //  * @brief Advance the electrolyte concentration by one timestep.
    //  *
    //  * Performs Crank–Nicolson assembly (left/right operators), applies
    //  * updated reaction and Neumann BC contributions, and solves for the
    //  * new true-DoF concentration vector before distributing back to \p Cn.
    //  *
    //  * @param Rx  Reaction field (input; used to assemble Rxe).
    //  * @param Cn  Electrolyte concentration grid function (in/out).
    //  * @param psx Phase field ψ used to mask diffusivity and BCs.
    //  */
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);


    // /**
    //  * @brief Advance the electrolyte concentration by one timestep.
    //  *
    //  * Performs Crank–Nicolson assembly (left/right operators), applies
    //  * updated reaction and Neumann BC contributions, and solves for the
    //  * new true-DoF concentration vector before distributing back to \p Cn.
    //  *
    //  * @param RxC  Reaction field cathode (input; used to assemble Rxe).
    //  * @param RxA  Reaction field anode (input; used to assemble Rxe).
    //  * @param Cn  Electrolyte concentration grid function (in/out).
    //  * @param psx Phase field ψ used to mask diffusivity and BCs.
    //  */
    void UpdateConcentration(mfem::ParGridFunction &RxC, mfem::ParGridFunction &RxA, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    // ---------- Geometry / spaces ----------
    Initialize_Geometry &geometry;                              ///< Geometry/mesh owner
    Domain_Parameters   &domain_parameters;                     ///< Domain/physics parameters
    BoundaryConditions   &boundary_conditions;
    sim::CellMode mode_;

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;       ///< Parallel FESpace
    FEMOperators fem;
    Utils utils;

    

    mfem::HypreParVector CeVn; ///< Concentration at next time step

    // ---------- Boundary forcing ----------
    mfem::ParLinearForm Fet; ///< Linear form for domain/boundary forcing
    mfem::Array<int> nbc_bdr; ///< Boundary markers for Neumann BCs
    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; ///< Product coefficient for boundary terms
    mfem::ConstantCoefficient nbcCoef; ///< Scalar coefficient in boundary term
    mfem::GridFunctionCoefficient matCoef_R; ///< Coefficient wrapper for PeR
    mfem::Array<int> boundary_dofs; ///< Boundary degrees of freedom


private:

    // ---------- Operators (forms and matrices) ----------
    std::unique_ptr<mfem::ParBilinearForm> Me_init; ///< Mass form
    std::unique_ptr<mfem::ParLinearForm>   Be_init; ///< Domain+boundary forcing form

    mfem::HypreParMatrix Mmate; ///< Mass matrix M
    mfem::HypreParMatrix Kmate; ///< Diffusion stiffness matrix K
    std::unique_ptr<mfem::ParBilinearForm> Ke2; ///< Stiffness form for diffusion

    // ---------- Solver / preconditioner ----------
    mfem::CGSolver      Me_solver; ///< Conjugate gradient solver
    mfem::HypreSmoother Me_prec;   ///< Hypre smoother preconditioner

    // ---------- Coefficients / fields ----------
    mfem::GridFunctionCoefficient cAe; ///< Reaction coefficient wrapper
    mfem::GridFunctionCoefficient cDe; ///< Diffusivity coefficient wrapper

    mfem::ParGridFunction De;  ///< Electrolyte diffusivity field
    mfem::ParGridFunction Rxe; ///< Reaction field
    mfem::ParGridFunction PeR; ///< Reaction potential field

    // ---------- State values ----------
    double eCrnt = 0.0; ///< Reaction current
    double infx  = 0.0; ///< Boundary flux value
    double L_w   = 0.0; ///< Characteristic electrolyte width

    // ---------- Vectors ----------
    mfem::HypreParVector Feb; ///< RHS vector
    mfem::HypreParVector X1v; ///< Scratch vector

    std::unique_ptr<mfem::HypreParMatrix> TmatR; ///< CN right operator (M - 0.5 dt K)
    std::unique_ptr<mfem::HypreParMatrix> TmatL; ///< CN left operator (M + 0.5 dt K)

    mfem::HypreParVector CeV0; ///< Current concentration (true DoFs)
    mfem::HypreParVector RHCe; ///< RHS after applying TmatR + forcingt


};

/** @} */ // end of group concentrations


#endif // CNE_HPP
