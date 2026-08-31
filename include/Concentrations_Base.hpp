#pragma once

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "BoundaryConditions.hpp"
#include "Utils.hpp"
#include "FEMOperators.hpp"
#include <memory>
#include <vector>

/**
 * @class ConcentrationBase
 * @brief Abstract base class for BESFEM concentration solvers.
 *
 * ConcentrationBase defines the shared interface and storage used by all
 * concentration-update modules, including solid-particle concentration and
 * electrolyte concentration solvers.
 *
 * Derived classes must implement:
 * - SetupField()
 * - UpdateConcentration()
 *
 * This class stores common mesh data, finite element spaces, element volumes,
 * diagnostic quantities, and temporary MFEM workspaces used during
 * concentration updates.
 */
class ConcentrationBase {
protected:
    /// Material type associated with this concentration solver.
    sim::MaterialType material;

    double net_pair_source = 0.0;
    double absolute_pair_source = 0.0;

public:
    /**
     * @brief Construct the base concentration solver.
     *
     * Stores references to geometry, domain parameters, and configuration data.
     * Also initializes shared mesh and finite-element-space data used by derived
     * concentration solvers.
     *
     * @param geo Geometry handler containing mesh and FE-space information.
     * @param para Domain parameters containing fields, volumes, and totals.
     * @param mat Material type associated with this solver.
     * @param cfg Simulation configuration settings.
     */
    ConcentrationBase(Initialize_Geometry &geo,
                      Domain_Parameters &para,
                      sim::MaterialType mat,
                      const SimulationConfig &cfg);

    /// Virtual destructor for safe cleanup through base-class pointers.
    virtual ~ConcentrationBase() = default;

    /**
     * @brief Initialize a concentration field and assemble model-specific data.
     *
     * Derived classes use this function to initialize the concentration field,
     * apply phase-field masks, and assemble any operators or coefficients needed
     * for the concentration update.
     *
     * @param Cn Concentration field to initialize.
     * @param initial_value Scalar initial concentration value.
     * @param psx Phase-field mask associated with this concentration field.
     * @param gtPsx Global normalization/integral of `psx`.
     */
    virtual void SetupField(mfem::ParGridFunction &Cn,
                            double initial_value,
                            mfem::ParGridFunction &psx,
                            double gtPsx) = 0;

    /**
     * @struct PairCoupling
     * @brief Workspace fields used for pairwise particle coupling.
     *
     * PairCoupling stores temporary fields needed when concentration updates
     * depend on interactions between neighboring particles, such as interfacial
     * coupling or chemical-potential exchange terms.
     */
    struct PairCoupling
    {
        mfem::ParGridFunction *sum_part = nullptr; ///< Accumulated pair contribution.
        mfem::ParGridFunction *weight   = nullptr; ///< Pair coupling weight field.
        mfem::ParGridFunction *grad_psi = nullptr; ///< Interface/gradient field between particles.
        mfem::ParGridFunction *mu_self  = nullptr; ///< Chemical potential of the current particle.
        mfem::ParGridFunction *mu_nbr   = nullptr; ///< Chemical potential of the neighboring particle.

        // int self_index = -1;
        // int neighbor_index = -1;
    };

    /**
     * @brief Advance the concentration field by one timestep.
     *
     * Derived classes implement the model-specific concentration update,
     * including reaction terms, diffusion or Cahn--Hilliard terms, masking, and
     * optional pairwise particle coupling.
     *
     * @param Rx Reaction/source term field.
     * @param Cn Concentration field to update.
     * @param psx Phase-field mask for the active region.
     * @param gtPsx Global normalization/integral of `psx`.
     * @param weight_elec Electrolyte/solid interface weighting field.
     * @param pair_terms Optional pairwise coupling workspaces.
     */
    virtual void UpdateConcentration(
        mfem::ParGridFunction &Rx,
        mfem::ParGridFunction &Cn,
        mfem::ParGridFunction &psx,
        double gtPsx,
        mfem::ParGridFunction &weight_elec,
        const std::vector<PairCoupling> &pair_terms) = 0;

    /**
     * @brief Apply salt or concentration conservation correction.
     *
     * The default implementation does nothing. Derived classes may override
     * this function when a conservation correction is needed after the
     * concentration update.
     *
     * @param Cn Concentration field to correct.
     * @param psx Phase-field mask defining the active region.
     */
    virtual void SaltConservation(mfem::ParGridFunction &Cn,
                                  mfem::ParGridFunction &psx)
    {}

    /**
     * @brief Return the current global lithiation or concentration fraction.
     *
     * The interpretation depends on the derived solver. For solid electrodes,
     * this is typically the average lithiation. For electrolyte solvers, it may
     * represent the normalized salt inventory.
     *
     * @return Global lithiation or concentration fraction.
     */
    double GetLithiation() const { return Xfr; }

    
    
    // double GetNetPairSource() const
    // {
    //     return net_pair_source;
    // }

    // double GetAbsolutePairSource() const
    // {
    //     return absolute_pair_source;
    // }


    // -------------------------------------------------------------------------
    // Global state / diagnostics
    // -------------------------------------------------------------------------
    double Xfr  = 0.0; ///< Global degree of lithiation or concentration.
    double CeC  = 0.0; ///< Local (per-element) salt or concentration moment.
    double gCeC = 0.0; ///< Global salt inventory (summed over domain).

    double CeAvg = 0.0; ///< Global average concentration.
    double Ce0   = 0.001; ///< Initial background concentration (reference).

    // -------------------------------------------------------------------------
    // Mesh / element / vertex information
    // -------------------------------------------------------------------------
    mfem::Array<double> VtxVal; ///< Vertex-value workspace buffer.

    int nE = 0; ///< Number of mesh elements.
    int nC = 0; ///< Number of nodes per element (corners).
    int nV = 0; ///< Total number of global vertices.

    mfem::Vector EAvg;          ///< Per-element averages (workspace).
    const mfem::Vector &EVol;   ///< Element volumes (from Domain_Parameters).

    // -------------------------------------------------------------------------
    // Phase-field normalization
    // -------------------------------------------------------------------------
    double gtPsi = 0.0; ///< Global ψ normalization (total).
    double gtPse = 0.0; ///< Global ψ normalization (electrolyte or solid).

    // -------------------------------------------------------------------------
    // Boundary and workspace DoF storage
    // -------------------------------------------------------------------------
    mfem::Array<int> boundary_dofs; ///< Boundary true DOFs.
    mfem::HypreParVector X1v;       ///< Scratch vector for assembly/solves.

    /// Reference to geometry data, including mesh and FE-space information.
    Initialize_Geometry &geometry;

    /// Reference to global simulation fields and domain parameters.
    Domain_Parameters &domain_parameters;

    /// Reference to user-defined simulation settings.
    const SimulationConfig &cfg;

    /// Shared parallel finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;

    /// Shared finite element helper operators.
    FEMOperators fem;

    /// Utility helper object for common postprocessing and calculations.
    Utils utils;

    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh (distributed).
    mfem::Mesh    *gmesh = nullptr; ///< Global serial mesh (rank 0, if available).

    std::unique_ptr<mfem::ParGridFunction> TmpF; ///< Temporary grid-function buffer.
};
