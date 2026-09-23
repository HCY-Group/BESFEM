#ifndef ELECTRODE_POTENTIAL_HPP
#define ELECTRODE_POTENTIAL_HPP

#include "Utils.hpp"
#include "Potentials_Base.hpp"
#include "SimTypes.hpp"

/**
 * @file ElectrodePotential.hpp
 * @brief Defines the solid-phase electric potential solver for electrode materials.
 */

class Initialize_Geometry;
class Domain_Parameters;
class BoundaryConditions;


/**
 * @class ElectrodePotential
 * @brief Solves the solid-phase electric potential within an electrode.
 *
 * ElectrodePotential assembles and solves the finite-element equations for
 * the solid-phase electric potential in either the anode or cathode. The
 * solver computes conductivity fields, incorporates Butler-Volmer reaction
 * source terms, applies boundary conditions, and updates the electrode
 * potential during each timestep.
 *
 * The class supports both single-particle and multi-particle electrode
 * simulations.
 */
class ElectrodePotential : public PotentialBase
{
public:
    
    /**
     * @brief Construct an electrode potential solver.
     *
     * @param geo Reference to the initialized geometry object.
     * @param para Reference to the domain parameter object.
     * @param bc Reference to the boundary-condition handler.
     * @param electrode Electrode being solved (anode or cathode).
     * @param material Active material associated with this electrode.
     * @param cfg Reference to the simulation configuration.
     */
    ElectrodePotential(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::Electrode electrode, sim::MaterialType material, const SimulationConfig &cfg);

    /**
     * @brief Initialize the electrode potential field.
     *
     * Sets the initial potential values and applies the electrode phase-field
     * mask prior to assembling the governing equations.
     *
     * @param ph Potential field to initialize.
     * @param initial_value Initial electrode potential.
     * @param psx Electrode phase-field mask.
     */
    void SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx);
    
    /**
     * @brief Solve the solid-phase potential equation.
     *
     * Assembles and solves the finite-element system for the electrode potential
     * using the current reaction source term and conductivity field.
     *
     * @param Rx Reaction/source term.
     * @param phx Electrode potential field.
     * @param psx Electrode phase-field mask.
     * @param gerror Output global solver error.
     */
    void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);

    /**
     * @brief Assemble the electrode potential system.
     *
     * Builds the conductivity matrix and reaction source terms for a single
     * electrode material.
     *
     * @param Cn Concentration field.
     * @param psx Electrode phase-field mask.
     * @param potential Electrode potential field.
     */
    void AssembleSystem(mfem::ParGridFunction &Cn,mfem::ParGridFunction &psx,mfem::ParGridFunction &potential) override;
    
    /**
     * @brief Assemble the electrode potential system for multiple particle groups.
     *
     * Computes conductivity contributions from multiple particle groups and
     * assembles the resulting finite-element system.
     *
     * @param Cn_groups Particle concentration fields.
     * @param psi_groups Particle phase-field masks.
     * @param materials Material associated with each particle group.
     * @param potential Electrode potential field.
     */
    void AssembleSystem(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups, const std::vector<sim::MaterialType> &materials, mfem::ParGridFunction &potential);

    /// Return the applied electrode boundary voltage.
    double GetBoundaryVoltage() const override { return Bv; }
    
    /// Increment the applied boundary voltage.
    void AddBoundaryVoltage(double dV){Bv += dV;}

    /// Return the GridFunction for conductivity
    mfem::ParGridFunction GetConductivity() override {return kap;}

private:

    Initialize_Geometry &geometry;          ///< Geometry and finite-element infrastructure.
    Domain_Parameters &domain_parameters;   ///< Domain fields and global quantities.
    BoundaryConditions &boundary_conditions;///< Boundary-condition manager.

    const SimulationConfig &cfg;            ///< Simulation settings.

    FEMOperators fem;                       ///< Finite-element assembly utilities.
    Utils utils;                            ///< Utility/helper functions.

    sim::Electrode electrode;               ///< Electrode being solved.
    sim::MaterialType material;             ///< Active material for this solver.

    mfem::ParGridFunction kap;              ///< Electrical conductivity field.
    mfem::GridFunctionCoefficient cKp;      ///< Conductivity coefficient.
    mfem::ParGridFunction RpP;              ///< Reaction source field.
    mfem::GridFunctionCoefficient cRp;      ///< Reaction coefficient.

    std::unique_ptr<mfem::ParLinearForm> Bp2; ///< Reaction linear form.
    std::unique_ptr<mfem::ParBilinearForm> Kp2; ///< Conductivity bilinear form.

    mfem::ParGridFunction pP0;              ///< Previous potential field.

    mfem::ParLinearForm B1t;                ///< Temporary linear form.
    mfem::HypreParVector X1v;               ///< Solution vector.
    mfem::HypreParVector B1v;               ///< Right-hand-side vector.
    mfem::HypreParVector Fpb;               ///< Reaction contribution vector.
    mfem::HypreParVector Xs0;               ///< Previous solution true-DOF vector.
    mfem::ParLinearForm Fpt;                ///< Reaction linear form.
    mfem::HypreParMatrix KmP;               ///< Potential system matrix.

    std::unique_ptr<mfem::HypreBoomerAMG> Mpp; ///< AMG preconditioner.

    mfem::Array<int> ess_tdof_list;         ///< Essential true DOFs.

    double Bv = 0.0;                        ///< Applied boundary voltage.
    double gtPs = 0.0;                      ///< Global phase-field normalization.

    mfem::CGSolver cgPP_solver;             ///< Conjugate-gradient solver.

    /**
     * @brief Compute the conductivity field for multiple particle groups.
     *
     * Combines the conductivity contributions from all particle groups into a
     * single effective conductivity field used during system assembly.
     *
     * @param Cn_groups Particle concentration fields.
     * @param psi_groups Particle phase-field masks.
     * @param materials Material associated with each particle group.
     */
    void ParticleConductivityMulti(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups, const std::vector<sim::MaterialType> &materials);
};

#endif // ELECTRODE_POTENTIAL_HPP
