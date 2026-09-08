#pragma once

/**
 * @file SimTypes.hpp
 * @brief Defines small enumeration types used throughout BESFEM simulations.
 *
 * Contains standard enumerations for cell configuration (half/full cell) and 
 * electrode selection (anode, cathode, or both in full-cell mode).
 */

namespace sim {

/**
 * @enum CellMode
 * @brief Specifies whether the simulation is a half-cell or full-cell.
 *
 * - **HALF** — A half-cell simulation with a single solid electrode 
 *              (anode or cathode) coupled to an electrolyte.
 * - **FULL** — A full-cell simulation including *both* electrodes 
 *              (anode + cathode) separated by an electrolyte domain.
 */
enum class CellMode {
    HALF,  ///< Half-cell configuration
    FULL   ///< Full-cell configuration
};

/**
 * @enum Electrode
 * @brief Identifies which electrode is active in a simulation.
 *
 * - **ANODE**   — Anode domain only  
 * - **CATHODE** — Cathode domain only  
 * - **BOTH**    — Used only for full-cell mode (anode + cathode)
 */
enum class Electrode {
    ANODE,    ///< Solid anode electrode
    CATHODE,  ///< Solid cathode electrode
    BOTH      ///< Both electrodes (valid only for full-cell simulations)
};


/**
 * @enum MaterialType
 * @brief Enumerates the types of electrode materials supported in the simulation.
 * 
 * - **Graphite** — Common anode material
 * - **NMC**       — Nickel Manganese Cobalt Oxide, a common cathode material
 * - **LFP**       — Lithium Iron Phosphate, another common cathode material
 * - **Electrolyte** — Represents the electrolyte domain (not a solid electrode)
 */

enum class MaterialType {
    Graphite, ///< Anode material
    NMC,       ///< Cathode material: Nickel Manganese Cobalt Oxide
    LFP,        ///< Cathode material: Lithium Iron Phosphate
    Carbon,   ///< Anode material: Carbon
    Electrolyte   ///< Electrolyte material
};

/**
 * @enum StopMode
 * @brief Defines the stopping condition for the simulation.
 * 
 * - **STEPS**   — Stop after a fixed number of timesteps
 * - **VOLTAGE** — Stop when the cell voltage reaches a specified threshold
 */
enum class StopMode
{
    STEPS,
    VOLTAGE
};


/**
 * @enum GeometryPhase
 * @brief Represents the phase of the geometry in the simulation.
 *
 * - **SOLID**     — Solid electrode phase
 * - **ELECTROLYTE** — Electrolyte phase
 */
enum class GeometryPhase
{
    SOLID,
    ELECTROLYTE
};

/**
 * @enum BoundarySide
 * @brief Represents the sides of the simulation domain for boundary conditions.
 *
 * - **WEST**   — West boundary
 * - **EAST**   — East boundary
 * - **SOUTH**  — South boundary
 * - **NORTH**  — North boundary
 * - **BOTTOM** — Bottom boundary (for 3D simulations)
 * - **TOP**    — Top boundary (for 3D simulations)
 */
enum class BoundarySide
{
    WEST,
    EAST,
    SOUTH,
    NORTH,
    BOTTOM,
    TOP
};


} // namespace sim
