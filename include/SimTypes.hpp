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

} // namespace sim
