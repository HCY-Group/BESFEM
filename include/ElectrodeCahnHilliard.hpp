#ifndef ELECTRODE_CAHN_HILLIARD_HPP
#define ELECTRODE_CAHN_HILLIARD_HPP

#include "Utils.hpp"
#include "Concentrations_Base.hpp"

/**
 * @file ElectrodeCahnHilliard.hpp
 * @brief Defines the Cahn--Hilliard concentration solver for electrode particles.
 */

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class ElectrodeCahnHilliard
 * @brief Solves solid-phase electrode concentration using a Cahn-Hilliard model.
 *
 * ElectrodeCahnHilliard updates the lithium concentration inside electrode
 * particles using a Cahn-Hilliard formulation with reaction coupling at the
 * particle/electrolyte interface.
 *
 * The class assembles the MFEM operators needed for the concentration solve,
 * evaluates chemical potential and mobility fields, and advances the particle
 * concentration each timestep.
 */
class ElectrodeCahnHilliard : public ConcentrationBase
{
public:
    /**
     * @brief Construct a Cahn-Hilliard electrode concentration solver.
     *
     * @param geo Reference to the initialized geometry object.
     * @param para Reference to the domain parameter object.
     * @param mat Material type for this electrode phase.
     * @param cfg Reference to the simulation configuration.
     */
    ElectrodeCahnHilliard(Initialize_Geometry &geo,
                          Domain_Parameters &para,
                          sim::MaterialType mat,
                          const SimulationConfig &cfg);

    /**
     * @brief Initialize the concentration field and assemble CH operators.
     *
     * Applies the phase-field mask, initializes the concentration field, and
     * assembles the mass, gradient, and auxiliary operators used by the
     * Cahn-Hilliard update.
     *
     * @param Cn Concentration field to initialize.
     * @param initial_value Initial lithium concentration.
     * @param psx Particle phase-field mask.
     * @param gtPsx Global integral of `psx`.
     */
    void SetupField(mfem::ParGridFunction &Cn,
                    double initial_value,
                    mfem::ParGridFunction &psx,
                    double gtPsx);

    /**
     * @brief Advance the electrode concentration by one timestep.
     *
     * Updates the concentration using the reaction source term, mobility field,
     * chemical-potential field, and optional particle-particle coupling terms.
     *
     * @param Rx Reaction/source term field.
     * @param Cn Concentration field to update.
     * @param psx Particle phase-field mask.
     * @param gtPsx Global integral of `psx`.
     * @param weight_elec Electrolyte-interface weighting field.
     * @param pair_terms Optional pairwise coupling terms between particles.
     */
    void UpdateConcentration(
        mfem::ParGridFunction &Rx,
        mfem::ParGridFunction &Cn,
        mfem::ParGridFunction &psx,
        double gtPsx,
        mfem::ParGridFunction &weight_elec,
        const std::vector<ConcentrationBase::PairCoupling> &pair_terms);

    /// Reference to user-defined simulation settings.
    const SimulationConfig &cfg;

    /// Return the GridFunction for diffusivity
    mfem::ParGridFunction GetDiffusivity() override {
        mfem::mfem_error("GetDiffusivity not implemented for CahnHilliard");
    }

private:
    mfem::ParGridFunction Mub; ///< Chemical-potential field.
    mfem::ParGridFunction Mob; ///< Mobility field.
    mfem::ParGridFunction RxA; ///< Reaction/source workspace field.

    mfem::ParMesh *pmesh; ///< Parallel mesh pointer.

    mfem::HypreParVector Lp1; ///< First auxiliary CH solve vector.
    mfem::HypreParVector Lp2; ///< Second auxiliary CH solve vector.
    mfem::HypreParVector MuV; ///< Chemical-potential true-DOF vector.

    mfem::HypreParVector PsVc; ///< Phase-field true-DOF vector.
    mfem::HypreParVector CpV0; ///< Previous concentration true-DOF vector.
    mfem::HypreParVector RHCp; ///< Right-hand side vector for concentration solve.
    mfem::HypreParVector CpVn; ///< Updated concentration true-DOF vector.

    std::unique_ptr<mfem::ParBilinearForm> M_init; ///< Mass-form assembly object.
    mfem::HypreParMatrix MmatCH; ///< Cahn--Hilliard mass matrix.

    std::unique_ptr<mfem::ParBilinearForm> Grad_MForm; ///< Mobility-gradient bilinear form.
    std::unique_ptr<mfem::ParBilinearForm> Grad_EForm; ///< Gradient-energy bilinear form.

    mfem::HypreParMatrix Grad_MM; ///< Mobility-weighted gradient matrix.
    mfem::HypreParMatrix Grad_EM; ///< Gradient-energy matrix.

    std::unique_ptr<mfem::ParLinearForm> B_init; ///< Linear-form assembly object.

    mfem::HypreParVector Fcb; ///< Boundary/source contribution vector.
    mfem::ParLinearForm Fct; ///< Linear form for source/reaction terms.

    mfem::GridFunctionCoefficient cAp; ///< Coefficient wrapper for reaction/source field.
    mfem::GridFunctionCoefficient cDp; ///< Coefficient wrapper for mobility field.

    std::unique_ptr<mfem::HypreParMatrix> TmatCH; ///< System matrix for the CH concentration update.

    double gtPsA = 0.0; ///< Global particle phase-field integral.
    double gtPsi = 0.0; ///< Global solid phase-field integral.

    mfem::CGSolver MCH_solver; ///< Linear solver for the CH system.
    mfem::HypreSmoother MCH_prec; ///< Preconditioner for the CH solver.

    bool combine_particle_groups = false; ///< Whether particle groups are solved as one combined group.
};

#endif // ELECTRODE_CAHN_HILLIARD_HPP
