#ifndef ELECTRODE_DIFFUSION_HPP
#define ELECTRODE_DIFFUSION_HPP

#include "Utils.hpp"
#include "Concentrations_Base.hpp"
#include "SimTypes.hpp"

/**
 * @file ElectrodeDiffusion.hpp
 * @brief Defines the diffusion-based concentration solver for electrode particles.
 */


class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class ElectrodeDiffusion
 * @brief Solves solid-state lithium diffusion within electrode particles.
 *
 * ElectrodeDiffusion updates the lithium concentration in electrode particles
 * using a diffusion equation with electrochemical reaction coupling at the
 * particle/electrolyte interface.
 *
 * The class assembles the finite-element operators required for the diffusion
 * solve, evaluates concentration-dependent material properties, and advances
 * the particle concentration at each timestep.
 */
class ElectrodeDiffusion : public ConcentrationBase
{
public:
    /**
     * @brief Construct a diffusion-based electrode concentration solver.
     *
     * @param geo Reference to the initialized geometry object.
     * @param para Reference to the domain parameter object.
     * @param mat Material type associated with this electrode.
     * @param cfg Reference to the simulation configuration.
     */
    ElectrodeDiffusion(Initialize_Geometry &geo, Domain_Parameters &para, sim::MaterialType mat, const SimulationConfig &cfg);

    /**
     * @brief Initialize the diffusion solver.
     *
     * Initializes the concentration field, applies the phase-field mask, and
     * assembles the finite-element operators required for the diffusion equation.
     *
     * @param Cn Concentration field to initialize.
     * @param initial_value Initial lithium concentration.
     * @param psx Particle phase-field mask.
     * @param gtPsx Global integral of the phase-field mask.
     */
    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx);
    
    /**
     * @brief Advance the concentration field by one timestep.
     *
     * Solves the diffusion equation, including electrochemical reaction source
     * terms and optional particle-particle coupling contributions.
     *
     * @param Rx Reaction/source term field.
     * @param Cn Concentration field to update.
     * @param psx Particle phase-field mask.
     * @param gtPsx Global integral of the phase-field mask.
     * @param weight_elec Electrolyte-interface weighting field.
     * @param pair_terms Pairwise particle coupling workspaces.
     */
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
    double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms);

    /// Reference to the simulation configuration.
    const SimulationConfig& cfg;

    /// Return the GridFunction for diffusivity
    mfem::ParGridFunction GetDiffusivity() override {return Dp;}

private:

    mfem::ParGridFunction Dp;   ///< Diffusivity field.
    mfem::ParGridFunction Rxn;  ///< Reaction/source field.

    mfem::HypreParVector PsVc;  ///< Phase-field true-DOF vector.
    mfem::HypreParVector CpV0;  ///< Previous concentration vector.
    mfem::HypreParVector RHCp;  ///< Right-hand side vector.
    mfem::HypreParVector CpVn;  ///< Updated concentration vector.

    std::unique_ptr<mfem::ParBilinearForm> Mt; ///< Mass bilinear form.
    mfem::HypreParMatrix Mmatp;                ///< Mass matrix.

    mfem::CGSolver Mp_solver;                 ///< Diffusion linear solver.
    mfem::HypreSmoother Mp_prec;              ///< Solver preconditioner.

    std::unique_ptr<mfem::ParBilinearForm> Kc2; ///< Diffusion bilinear form.
    mfem::HypreParMatrix Kmatp;                 ///< Diffusion stiffness matrix.

    std::unique_ptr<mfem::ParLinearForm> Bc2; ///< Linear-form assembly.

    mfem::HypreParVector Fcb; ///< Boundary/source contribution vector.
    mfem::ParLinearForm Fct;  ///< Source-term linear form.

    mfem::GridFunctionCoefficient cAp; ///< Reaction-field coefficient.
    mfem::GridFunctionCoefficient cDp; ///< Diffusivity coefficient.

    std::unique_ptr<mfem::HypreParMatrix> Tmatp; ///< System matrix.

    double gtPsC = 0.0; ///< Global particle phase-field integral.
    double gtPsi = 0.0; ///< Global solid phase-field integral.

    bool combine_particle_groups = false; ///< Combine all particles into a single diffusion solve.

};

#endif // ELECTRODE_DIFFUSION_HPP
