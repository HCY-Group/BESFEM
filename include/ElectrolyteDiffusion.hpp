#ifndef ELECTROLYTE_DIFFUSION_HPP
#define ELECTROLYTE_DIFFUSION_HPP

#include "Utils.hpp"
#include "Concentrations_Base.hpp"

/**
 * @file ElectrolyteDiffusion.hpp
 * @brief Defines the electrolyte salt-concentration diffusion solver.
 */

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class ElectrolyteDiffusion
 * @brief Solves electrolyte salt transport in BESFEM battery simulations.
 *
 * ElectrolyteDiffusion updates the electrolyte concentration field using a
 * diffusion equation with electrochemical reaction coupling at electrode/
 * electrolyte interfaces. It also provides a salt-conservation correction used
 * to maintain the electrolyte inventory during time integration.
 */
class ElectrolyteDiffusion : public ConcentrationBase
{ 
public:

    /**
     * @brief Construct an electrolyte diffusion solver.
     *
     * @param geo Reference to the initialized geometry object.
     * @param para Reference to the domain parameter object.
     * @param bc Reference to the boundary-condition handler.
     * @param mode Cell mode, such as half-cell or full-cell.
     * @param mat Material type associated with this electrolyte solve.
     * @param cfg Reference to the simulation configuration.
     */
    ElectrolyteDiffusion(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::CellMode mode, sim::MaterialType mat, const SimulationConfig &cfg);

    /**
     * @brief Initialize the electrolyte concentration field and operators.
     *
     * Initializes the electrolyte concentration, applies the electrolyte phase-field
     * mask, and assembles the mass/diffusion operators used for the electrolyte
     * update.
     *
     * @param Cn Electrolyte concentration field to initialize.
     * @param initial_value Initial electrolyte concentration.
     * @param psx Electrolyte phase-field mask.
     * @param gtPsx Global integral of the electrolyte phase-field mask.
     */
    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx);
    
    /**
     * @brief Advance the electrolyte concentration by one timestep.
     *
     * Updates the electrolyte concentration using diffusion, reaction source terms,
     * boundary contributions, and the electrolyte phase-field mask.
     *
     * @param Rx Reaction/source term field.
     * @param Cn Electrolyte concentration field to update.
     * @param psx Electrolyte phase-field mask.
     * @param gtPsx Global integral of the electrolyte phase-field mask.
     * @param weight_elec Electrolyte-interface weighting field.
     * @param pair_terms Pairwise particle coupling terms; unused by the default electrolyte update.
     */
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                             double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms);

    /**
     * @brief Apply salt-conservation correction to the electrolyte concentration.
     *
     * Corrects the electrolyte concentration field so the total electrolyte salt
     * inventory remains consistent after the diffusion/reaction update.
     *
     * @param Cn Electrolyte concentration field to correct.
     * @param psx Electrolyte phase-field mask.
     */
    void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) override;

    /// Return the GridFunction for diffusivity
    mfem::ParGridFunction GetDiffusivity() override {return De;}

    BoundaryConditions &boundary_conditions; ///< Boundary-condition manager.
    sim::CellMode mode_;                     ///< Cell mode for this electrolyte solve.
    const SimulationConfig &cfg;             ///< Simulation configuration settings.

    mfem::HypreParVector CeVn;               ///< Updated electrolyte concentration true-DOF vector.
    mfem::ParLinearForm Fet;                 ///< Electrolyte boundary/source linear form.
    mfem::Array<int> nbc_bdr;                ///< Natural-boundary-condition marker array.
    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; ///< Natural-boundary-condition coefficient.
    mfem::ConstantCoefficient nbcCoef;       ///< Constant natural-boundary coefficient.
    mfem::GridFunctionCoefficient matCoef_R; ///< Reaction/source coefficient.
    mfem::Array<int> boundary_dofs;          ///< Boundary true DOFs.         


private:

    std::unique_ptr<mfem::ParBilinearForm> Me_init; ///< Electrolyte mass bilinear form.
    std::unique_ptr<mfem::ParLinearForm> Be_init;   ///< Electrolyte initial-condition linear form.

    mfem::HypreParMatrix Mmate; ///< Electrolyte mass matrix.
    mfem::HypreParMatrix Kmate; ///< Electrolyte diffusion stiffness matrix.
    std::unique_ptr<mfem::ParBilinearForm> Ke2; ///< Electrolyte diffusion bilinear form.

    mfem::CGSolver Me_solver;      ///< Electrolyte linear solver.
    mfem::HypreSmoother Me_prec;   ///< Electrolyte solver preconditioner.

    mfem::GridFunctionCoefficient cAe; ///< Electrolyte reaction/source coefficient.
    mfem::GridFunctionCoefficient cDe; ///< Electrolyte diffusivity coefficient.

    mfem::ParGridFunction De;  ///< Electrolyte diffusivity field.
    mfem::ParGridFunction Rxe; ///< Electrolyte reaction/source field.
    mfem::ParGridFunction PeR; ///< Electrolyte phase/reaction workspace field.

    double eCrnt = 0.0; ///< Electrolyte current/source diagnostic.
    double infx  = 0.0; ///< Boundary flux diagnostic.
    double L_w   = 0.0; ///< Salt-conservation correction factor.

    mfem::HypreParVector Feb; ///< Boundary/source contribution vector.
    mfem::HypreParVector X1v; ///< Temporary solution/work vector.

    std::unique_ptr<mfem::HypreParMatrix> TmatR; ///< Right-hand-side update matrix.
    std::unique_ptr<mfem::HypreParMatrix> TmatL; ///< Left-hand-side system matrix.

    mfem::HypreParVector CeV0; ///< Previous electrolyte concentration true-DOF vector.
    mfem::HypreParVector RHCe; ///< Electrolyte right-hand-side vector.
    mfem::HypreParVector PsVc; ///< Electrolyte phase-field true-DOF vector.

};

#endif // ELECTROLYTE_DIFFUSION_HPP
