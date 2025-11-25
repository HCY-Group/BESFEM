#ifndef REACTION_HPP
#define REACTION_HPP

#include "Concentrations_Base.hpp"

/**
 * @file Reaction.hpp
 * @brief Electrochemical reaction model for BESFEM.
 *
 * Implements Butler–Volmer kinetics, exchange current density evaluation,
 * open-circuit voltage (OCV) lookup, forward/backward reaction rate
 * constants, and global reaction current integration. Supports
 * single-interface (half-cell) and multi-interface (full-cell)
 * configurations.
 */

/**
 * @class Reaction
 * @brief Electrochemical reaction operator (Butler–Volmer kinetics).
 *
 * Responsibilities:
 * - Initialize reaction fields
 * - Compute exchange-current density from lookup tables
 * - Compute forward/backward rate constants
 * - Evaluate Butler–Volmer reaction rates on interfaces
 * - Compute global reaction current (MPI reduced)
 *
 * Supports both:
 * - Half-cell models (single interface)
 * - Full-cell models (anode + cathode + electrolyte coupling)
 */
class Reaction {
public:

    /**
     * @brief Construct the Reaction operator.
     *
     * Stores references to geometry and domain parameters and initializes
     * workspace fields.
     *
     * @param geo  Geometry/mesh handler.
     * @param para Domain parameter container.
     */
    Reaction(Initialize_Geometry &geo, Domain_Parameters &para);

    Initialize_Geometry &geometry;          ///< Geometry and mesh infrastructure.
    Domain_Parameters   &domain_parameters; ///< Material and phase-field parameters.

    /**
     * @brief Assign an initial reaction rate.
     *
     * @param Rx            Reaction-rate field (in/out).
     * @param initial_value Scalar initialization value.
     */
    void Initialize(mfem::ParGridFunction &Rx, double initial_value);

    /**
     * @brief Compute exchange-current density using lookup tables.
     *
     * Uses concentration-dependent \( i_0(c) \), OCV(c), and mobility tables
     * to populate reaction-related material fields.
     *
     * @param Cn Concentration field used in the lookup.
     */
    void TableExchangeCurrentDensity(mfem::ParGridFunction &Cn);

    /**
     * @brief Compute exchange-current density (single concentration field).
     *
     * @param Cn Concentration field.
     */
    void ExchangeCurrentDensity(mfem::ParGridFunction &Cn);

    /**
     * @brief Compute exchange-current density (two concentration fields).
     *
     * Used for full-cell reaction evaluation (e.g., anode/electrolyte interface).
     *
     * @param Cn1 Concentration on one interface side.
     * @param Cn2 Concentration on the other interface side.
     */
    void ExchangeCurrentDensity(mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2);

    /**
     * @brief Apply the Butler–Volmer equation for a single interface.
     *
     * Computes using concentration and potential fields.
     *
     * @param Rx   Output reaction-rate field (in/out).
     * @param Cn1  Concentration on interface side 1.
     * @param Cn2  Concentration on interface side 2.
     * @param phx1 Potential on interface side 1.
     * @param phx2 Potential on interface side 2.
     */
    void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1,
                      mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);

    /**
     * @brief Butler–Volmer update for full-cell (3-phase) systems.
     *
     * Computes reaction rates for:
     * - cathode/electrolyte interface
     * - anode/electrolyte interface
     * - the combined reaction field Rx
     *
     * @param Rx   Combined reaction-rate field (in/out).
     * @param Rx1  Cathode-side reaction field.
     * @param Rx2  Anode-side reaction field.
     * @param Cn1  Cathode-side concentration.
     * @param Cn2  Anode-side concentration.
     * @param Cn3  Electrolyte concentration.
     * @param phx1 Cathode potential.
     * @param phx2 Anode potential.
     * @param phx3 Electrolyte potential.
     */
    void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2,
                      mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &Cn3,
                      mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2, mfem::ParGridFunction &phx3);

    /**
     * @brief Compute global reaction current.
     *
     * Integrates \( R(x)\,A_v(x) \) across the domain using MPI reduction.
     *
     * @param Rx             Reaction-rate field.
     * @param global_current Output: MPI-reduced total reaction current.
     */
    void TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current);

    // -------------------------------------------------------------------------
    // Material and reaction-related grid functions (owned)
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParGridFunction> Kfw;  ///< Forward rate constant k_f.
    std::unique_ptr<mfem::ParGridFunction> Kbw;  ///< Backward rate constant k_b.
    std::unique_ptr<mfem::ParGridFunction> KfA;  ///< Forward rate (anode).
    std::unique_ptr<mfem::ParGridFunction> KbA;  ///< Backward rate (anode).
    std::unique_ptr<mfem::ParGridFunction> KfC;  ///< Forward rate (cathode).
    std::unique_ptr<mfem::ParGridFunction> KbC;  ///< Backward rate (cathode).

    std::unique_ptr<mfem::ParGridFunction> i0C;  ///< Exchange-current density.
    std::unique_ptr<mfem::ParGridFunction> OCV;  ///< Open-circuit voltage.

    std::unique_ptr<mfem::ParGridFunction> i0CA; ///< Exchange-current density (anode).
    std::unique_ptr<mfem::ParGridFunction> OCVA; ///< OCV table (anode).
    std::unique_ptr<mfem::ParGridFunction> i0CC; ///< Exchange-current density (cathode).
    std::unique_ptr<mfem::ParGridFunction> OCVC; ///< OCV table (cathode).

    // -------------------------------------------------------------------------
    // Non-owning surface-area fields
    // -------------------------------------------------------------------------
    const mfem::ParGridFunction *AvP = nullptr; ///< Particle surface area.
    const mfem::ParGridFunction *AvB = nullptr; ///< Boundary surface area.
    const mfem::ParGridFunction *AvA = nullptr; ///< Anode surface area.
    const mfem::ParGridFunction *AvC = nullptr; ///< Cathode surface area.

private:

    /**
     * @brief Set the initial reaction rate value.
     *
     * Internal helper called by Initialize().
     *
     * @param Cn            Concentration field.
     * @param initial_value Initial scalar value.
     */
    void SetInitialReaction(mfem::ParGridFunction &Cn, double initial_value);

    // -------------------------------------------------------------------------
    // Mesh and FE space
    // -------------------------------------------------------------------------
    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh (non-owning).
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

    // -------------------------------------------------------------------------
    // Mesh-related parameters
    // -------------------------------------------------------------------------
    int nE = 0; ///< Number of elements.
    int nC = 0; ///< Number of corners per element.
    int nV = 0; ///< Number of vertices.

    double local_current = 0.0; ///< Local reaction current (before MPI reduce).

    std::unique_ptr<mfem::ParGridFunction> dPHE; ///< Potential drop (electrolyte).
    std::unique_ptr<mfem::ParGridFunction> dPHA; ///< Potential drop (anode).
    std::unique_ptr<mfem::ParGridFunction> dPHC; ///< Potential drop (cathode).

    const mfem::Vector &EVol; ///< Element volumes.

    // -------------------------------------------------------------------------
    // Lookup-table data (size = 101)
    // -------------------------------------------------------------------------
    mfem::Vector Ticks     = mfem::Vector(101); ///< Concentration ticks.
    mfem::Vector chmPot    = mfem::Vector(101); ///< Charge-transfer potential.
    mfem::Vector Mobility  = mfem::Vector(101); ///< Mobility table.
    mfem::Vector OCV_file  = mfem::Vector(101); ///< OCV lookup.
    mfem::Vector i0_file   = mfem::Vector(101); ///< Exchange-current lookup.

    /**
     * @brief Linearly interpolate from a reaction lookup table.
     *
     * @param cn    Concentration value (0–1).
     * @param ticks Sorted tick positions.
     * @param data  Tabulated values.
     * @return Interpolated value.
     */
    double GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data);
};

#endif // REACTION_HPP
