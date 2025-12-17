#ifndef CNC_HPP
#define CNC_HPP

#include "Concentrations_Base.hpp"

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class CnC
 * @brief Diffusion-based concentration solver for cathode (or generic solid) particles.
 *
 * The CnC class implements a standard diffusion model
 * for solid active-material particles. It assembles mass and stiffness
 * operators, applies diffusivity and reaction coefficients, and updates the
 * concentration field using implicit/explicit time-stepping based on MFEM
 * and Hypre.
 */
class CnC : public ConcentrationBase {
public:

    /**
     * @brief Construct the diffusion-based concentration solver.
     *
     * @param geo  Geometry handler (parallel mesh, FE space, DOF mappings).
     * @param para Domain parameters (material properties, operating mode).
     */
    CnC(Initialize_Geometry &geo, Domain_Parameters &para);

    /**
     * @brief Initialize concentration field and assemble the diffusion operators.
     *
     * Performs:
     * - Initial concentration assignment (masked by ψ)
     * - Mass-matrix assembly
     * - Diffusion stiffness-matrix assembly
     * - Setup of CG solver and preconditioner
     *
     * @param Cn            Concentration field to initialize.
     * @param initial_value Initial scalar value for the concentration.
     * @param psx           Phase-field ψ identifying solid region.
     */
    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Advance concentration by one timestep using diffusion.
     *
     * Uses current reaction field, diffusivity, and stiffness/mass operators
     * to compute the updated concentration. Applies ψ-field masking to restrict
     * the update to the solid region.
     *
     * @param Rx  Reaction term applied to the concentration.
     * @param Cn  Concentration field (input/output).
     * @param psx Phase-field ψ for masking solid region.
     */
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

private:

    // -------------------------------------------------------------------------
    // Geometry / Spaces
    // -------------------------------------------------------------------------
    Initialize_Geometry &geometry;        ///< Geometry and mesh structure.
    Domain_Parameters   &domain_parameters; ///< Material and physics parameters.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

    FEMOperators fem;  ///< High-level operator builder (mass, stiffness, forms).
    Utils utils;       ///< Utility helpers (clamping, interpolation, etc.).

    // -------------------------------------------------------------------------
    // Grid-function state (continuous FEM fields)
    // -------------------------------------------------------------------------
    mfem::ParGridFunction RxC; ///< Reaction term mapped to the concentration space.
    mfem::ParGridFunction Dp;  ///< Diffusivity field D(c) over the solid region.

    // -------------------------------------------------------------------------
    // True DoF vectors (Hypre Parallel)
    // -------------------------------------------------------------------------
    mfem::HypreParVector PsVc; ///< ψ-field mask in true-DoF representation.
    mfem::HypreParVector CpV0; ///< Concentration vector at the current timestep.
    mfem::HypreParVector RHCp; ///< RHS contributions (mass-matrix product + reaction).
    mfem::HypreParVector CpVn; ///< Concentration vector at next timestep.

    // -------------------------------------------------------------------------
    // Operators (mass/stiffness/forcing)
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParBilinearForm> Mt; ///< Mass bilinear form.
    mfem::HypreParMatrix Mmatp;                ///< Mass matrix.

    mfem::CGSolver       Mp_solver;            ///< Solver for mass solves M x = b.
    mfem::HypreSmoother  Mp_prec;              ///< Preconditioner for Mp_solver.

    std::unique_ptr<mfem::ParBilinearForm> Kc2; ///< Diffusion stiffness bilinear form.
    mfem::HypreParMatrix Kmatp;                 ///< Diffusion stiffness matrix.

    std::unique_ptr<mfem::ParLinearForm> Bc2;   ///< Forcing linear form.

    mfem::HypreParVector Fcb; ///< RHS force-vector (true DoFs).
    mfem::ParLinearForm Fct; ///< Assembled linear form for forcing.

    // -------------------------------------------------------------------------
    // Coefficients (wrap grid functions into FE operators)
    // -------------------------------------------------------------------------
    mfem::GridFunctionCoefficient cAp; ///< Reaction coefficient wrapper.
    mfem::GridFunctionCoefficient cDp; ///< Diffusivity coefficient wrapper.

    // -------------------------------------------------------------------------
    // System Matrix for Time-Stepping
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::HypreParMatrix> Tmatp; ///< Implicit time-step system matrix.

    // -------------------------------------------------------------------------
    // Global Scaling
    // -------------------------------------------------------------------------
    double gtPsC = 0.0; ///< Global normalization for ψ (solid).
    double gtPsi = 0.0; ///< Global normalization for ψ (total).
};

#endif // CNC_HPP
