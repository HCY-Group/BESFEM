#ifndef CNA_HPP
#define CNA_HPP

#include "mfem.hpp"
#include "Concentrations_Base.hpp"
#include <memory>
#include <vector>

/**
 * @file CnA.hpp
 * @brief Cahn–Hilliard concentration evolution for the anode (solid phase).
 *
 * Defines the CnA class, which advances a solid-state concentration field
 * using a Cahn–Hilliard formulation, including mobility-weighted diffusion,
 * chemical potential calculation, reaction terms, and MFEM/Hypre operators.
 */

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class CnA
 * @brief Cahn–Hilliard concentration solver for the solid electrode phase.
 *
 * Implements a fully-coupled Cahn–Hilliard update for the concentration field
 * inside the solid active material of the anode.  
 *
 */
class CnA : public ConcentrationBase {
public:

    /**
     * @brief Construct the Cahn–Hilliard solver for the solid anode.
     *
     * @param geo  Geometry handler (mesh, FE space, DOF mappings).
     * @param para Domain and physics parameters.
     */
    CnA(Initialize_Geometry &geo, Domain_Parameters &para);

    /**
     * @brief Initialize concentration state and assemble operators.
     *
     * Performs:
     * - Initial concentration assignment.
     * - Initial Lithiation fraction computation (masked by ψ).
     * - Mass/stiffness matrix assembly.
     * - CG solver + Hypre preconditioner setup.
     *
     * @param Cn            Concentration grid function to initialize.
     * @param initial_value Scalar initial value for the concentration field.
     * @param psx           Phase-field mask ψ defining the solid region.
     */
    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Advance concentration by one timestep using the Cahn–Hilliard update.
     *
     * Recomputes chemical potential μ(c), mobility M(c), and reaction contributions.
     * Applies masking via ψ and clamps concentration in void regions.
     *
     * @param Rx  Reaction term (mapped to each DoF).
     * @param Cn  Concentration field (input/output).
     * @param psx Phase-field mask ψ for solid-region weighting.
     */
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

private:
    // -------------------------------------------------------------------------
    // Geometry and spaces
    // -------------------------------------------------------------------------
    Initialize_Geometry &geometry;     ///< Geometry/mesh container.
    Domain_Parameters   &domain_parameters; ///< Physics and material parameters.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.
    mfem::ParMesh *pmesh;              ///< Parallel mesh pointer.

    FEMOperators fem;                  ///< Matrix assembly helpers.
    Utils utils;                       ///< Utility helpers (clamping, etc.).

    // -------------------------------------------------------------------------
    // Grid-function state fields (continuous FEM)
    // -------------------------------------------------------------------------
    mfem::ParGridFunction Mub;   ///< Chemical potential μ(c).
    mfem::ParGridFunction Mob;   ///< Mobility field M(c)·ψ.
    mfem::ParGridFunction RxA;   ///< Reaction source mapped to the grid.

    // -------------------------------------------------------------------------
    // True-DoF vectors (Hypre parallel vectors)
    // -------------------------------------------------------------------------
    mfem::HypreParVector Lp1;    ///< Laplacian(C) intermediate.
    mfem::HypreParVector Lp2;    ///< Accumulated RHS contributions.
    mfem::HypreParVector MuV;    ///< μ(c) in true-DoF representation.
    mfem::HypreParVector PsVc;   ///< Phase-field ψ in true-DoF representation.

    mfem::HypreParVector CpV0;   ///< Concentration (current step).
    mfem::HypreParVector RHCp;   ///< Mass-matrix product + reaction RHS.
    mfem::HypreParVector CpVn;   ///< Previous-step concentration (optional).

    // -------------------------------------------------------------------------
    // Bilinear and linear forms (operators)
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParBilinearForm> M_init; ///< Mass form (assembled once).
    mfem::HypreParMatrix MmatCH;                    ///< Mass matrix M.

    mfem::CGSolver       MCH_solver;                ///< CG solver for M x = b.
    mfem::HypreSmoother  MCH_prec;                  ///< Preconditioner (Jacobi).

    std::unique_ptr<mfem::ParBilinearForm> Grad_MForm; ///< Mobility-weighted stiffness form.
    mfem::HypreParMatrix Grad_MM;                      ///< Mobility-weighted stiffness matrix.

    std::unique_ptr<mfem::ParLinearForm> B_init;       ///< Reaction linear form.

    std::unique_ptr<mfem::ParBilinearForm> Grad_EForm; ///< Energy Laplacian form.
    mfem::HypreParMatrix Grad_EM;                      ///< Energy operator matrix.

    mfem::ParLinearForm Fct;   ///< Free-energy / reaction contributions.
    mfem::HypreParVector Fcb;  ///< Assembled linear-form values (true DoFs).

    // -------------------------------------------------------------------------
    // Coefficients
    // -------------------------------------------------------------------------
    mfem::GridFunctionCoefficient cDp; ///< Mobility coefficient from Mob.
    mfem::GridFunctionCoefficient cAp; ///< Reaction coefficient from RxA.

    // -------------------------------------------------------------------------
    // Tabulated material data (size = 101)
    // -------------------------------------------------------------------------
    mfem::Vector Ticks   = mfem::Vector(101); ///< Concentration ticks (0.00 → 1.00).
    mfem::Vector chmPot  = mfem::Vector(101); ///< Chemical potential μ(c) table.
    mfem::Vector Mobility= mfem::Vector(101); ///< Mobility M(c) table.
    mfem::Vector OCV     = mfem::Vector(101); ///< Open-circuit voltage table.
    mfem::Vector i0      = mfem::Vector(101); ///< Exchange-current-density table.

    double gtPsA = 0.0; ///< Global normalization for ψ (solid-region weight).
    double gtPsi = 0.0; ///< Global normalization for ψ (total).

    /**
     * @brief Linearly interpolate tabulated material data at concentration `cn`.
     *
     * Clamps `cn` to (1e-6, 0.999999) and performs linear interpolation over
     * uniformly spaced concentration ticks (0.00, 0.01, ..., 1.00).
     *
     * @param cn    Concentration value in [0,1].
     * @param ticks Tick vector (monotonic, uniform spacing).
     * @param data  Data table aligned with `ticks`.
     * @return Interpolated value at concentration `cn`.
     */
    double GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data);
};

#endif // CNA_HPP
