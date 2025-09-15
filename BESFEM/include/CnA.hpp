// CnA.hpp
#ifndef CNA_HPP
#define CNA_HPP

#include "mfem.hpp"
#include "Concentrations_Base.hpp"
#include <memory>
#include <vector>

/**
 * @file CnA.hpp
 * @brief Cahn–Hilliard concentration evolution for battery simulations.
 *
 * This class advances a concentration field with a Cahn–Hilliard formulation.
 */

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @defgroup concentrations Concentration Modules
 * @brief Classes that advance concentration fields in the simulation.
 * @{
 */

/**
 * @class CnA
 * @brief Cahn–Hilliard concentration solver (solid electrode phase).
 * @ingroup concentrations
 *
 * @details
 * Implements a Cahn–Hilliard update for the concentration field:
 * - Builds/maintains mass and stiffness operators in parallel (MFEM/Hypre).
 * - Uses grid-function coefficients for mobility and reaction terms.
 * - Pulls chemical potential and mobility from tabulated data and linearly
 *   interpolates them per-DoF.
 *
 * Typical usage:
 * 1. Construct with geometry and domain parameters.
 * 2. Call Initialize() with the initial concentration and phase field \p psx.
 * 3. Repeatedly call TimeStep() each timestep.
 */

class CnA : public Concentrations {
public:
    /**
     * @brief Construct the Cahn–Hilliard concentration solver.
     *
     * @param geo  Geometry/space container (mesh, parallel FESpace, boundary DOFs).
     * @param para Domain/physics parameters (time step, material scales, etc.).
     */
    CnA(Initialize_Geometry &geo, Domain_Parameters &para);

    /**
     * @brief Assemble operators and prepare the solver.
     *
     * Sets the initial concentration value, computes initial lithiation from \p psx,
     * assembles the mass matrix and stiffness operators, and configures the CG
     * solver and preconditioner.
     *
     * @param Cn            Concentration field (solid) to initialize.
     * @param initial_value Initial scalar value used for Cn before masking.
     * @param psx           Phase field ψ (used for masking/region weighting).
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Advance the concentration by one timestep (Cahn–Hilliard update).
     *
     * Uses tabulated \e μ(C) and \e M(C) (mobility) to build the linear systems,
     * solves for the new concentration, and clamps void-region values.
     *
     * @param Rx  Reaction source field (input field used to assemble \e Rxc).
     * @param Cn  Concentration field (in/out).
     * @param psx Phase field ψ (masking/solid-region weighting).
     */
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    
private:
    // ---------- Geometry / spaces ----------
    Initialize_Geometry &geometry;                              ///< Geometry/mesh owner
    Domain_Parameters   &domain_parameters;                     ///< Domain/physics parameters
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;       ///< Parallel FESpace 
    mfem::ParMesh *pmesh;                                       ///< Parallel mesh

    // ---------- State fields (grid functions) ----------
    mfem::ParGridFunction Mub;   ///< Chemical potential field μ(C) 
    mfem::ParGridFunction Mob;   ///< Mobility field M(C)·ψ 
    mfem::ParGridFunction Rxc;   ///< Reaction source term mapped to grid

    // ---------- Work vectors (true DoFs / Hypre) ----------
    mfem::HypreParVector phV0;   ///< Generic scratch vector.
    mfem::HypreParVector Lp1;    ///< Laplacian( C ) intermediate.
    mfem::HypreParVector Lp2;    ///< Accumulated RHS contributions.
    mfem::HypreParVector RHS;    ///< Generic RHS buffer.
    mfem::HypreParVector MuV;    ///< True-DoF chemical potential vector.
    mfem::HypreParVector PsVc;   ///< True-DoF phase-field mask ψ.

    mfem::HypreParVector CpV0;   ///< Current concentration (true DoFs).
    mfem::HypreParVector RHCp;   ///< Mass-matrix product + source (RHS).
    mfem::HypreParVector CpVn;   ///< Previous-step concentration (optional).

    // ---------- Operators (forms and matrices) ----------
    std::unique_ptr<mfem::ParBilinearForm> M_init;  ///< Mass form (assembled once, updated if needed).
    mfem::HypreParMatrix MmatCH;                    ///< Mass matrix M for Cahn–Hilliard.

    mfem::CGSolver       MCH_solver;                ///< CG solver for the mass solve M x = b.
    mfem::HypreSmoother  MCH_prec;                  ///< Jacobi smoother (or other) as preconditioner.

    std::unique_ptr<mfem::ParBilinearForm> Grad_MForm; ///< Stiffness form for mobility-weighted operator.
    mfem::HypreParMatrix Grad_MM;                      ///< Stiffness matrix with mobility (K_M).

    std::unique_ptr<mfem::ParLinearForm>   B_init;     ///< Linear form for reaction/forcing terms.

    std::unique_ptr<mfem::ParBilinearForm> Grad_EForm; ///< Stiffness form for energy operator.
    mfem::HypreParMatrix Grad_EM;                      ///< Energy matrix (e.g., Laplacian scaling).

    mfem::ParLinearForm Fct;  ///< Assembled linear form for free-energy / reaction contributions.
    mfem::HypreParVector Fcb; ///< True-DoF vector holding assembled linear-form values.

    // ---------- Coefficients ----------
    mfem::GridFunctionCoefficient cDp; ///< Mobility coefficient wrapper (binds to @ref Mob).
    mfem::GridFunctionCoefficient cAp; ///< Reaction coefficient wrapper (binds to @ref Rxc).

    // ---------- Tabulated material data (size = 101) ----------
    mfem::Vector Ticks   = mfem::Vector(101); ///< Concentration ticks in [0,1].
    mfem::Vector chmPot  = mfem::Vector(101); ///< Chemical potential μ(C) table.
    mfem::Vector Mobility= mfem::Vector(101); ///< Mobility M(C) table (scaled).
    mfem::Vector OCV     = mfem::Vector(101); ///< Open-circuit voltage table.
    mfem::Vector i0      = mfem::Vector(101); ///< Exchange-current-density table.

    /**
     * @brief Linear interpolation utility for tabulated data.
     *
     * Bounds \p cn to (1e-6, 0.999999) and performs linear interpolation on
     * the provided \p data using uniform tick spacing (0.01).
     *
     * @param cn    Concentration value in [0,1].
     * @param ticks Tick vector (expected 0.00, 0.01, ..., 1.00).
     * @param data  Data values aligned with \p ticks.
     * @return Interpolated value at concentration \p cn.
     */    
    double GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data);


};

/** @} */ // end of group concentrations


#endif // CNA_HPP
