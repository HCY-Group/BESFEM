#pragma once

/**
 * @file SimTypes.hpp
 * @brief Defines small enumeration types used throughout BESFEM simulations.
 *
 * Contains standard enumerations for cell configuration (half/full cell) and 
 * electrode selection (anode, cathode, or both in full-cell mode).
 */

/**
 * @namespace sim
 * @brief Contains simulation types used throughout BESFEM.
 */
namespace sim {

/**
 * @enum CellMode
 * @brief Specifies whether the simulation is a half-cell or full-cell.
 *
 */
enum class CellMode {
    HALF,  ///< Half-cell configuration
    FULL   ///< Full-cell configuration
};

/**
 * @enum Electrode
 * @brief Identifies which electrode is active in a simulation.
 * 
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
 */

enum class MaterialType {
    Graphite, ///< Anode material
    NMC,       ///< Cathode material: Nickel Manganese Cobalt Oxide
    LFP,        ///< Cathode material: Lithium Iron Phosphate
    Carbon,   ///< Anode material: Carbon
    Silicon,  ///< Anode material: Silicon
    Electrolyte   ///< Electrolyte material
};

/**
 * @enum StopMode
 * @brief Defines the stopping condition for the simulation.
 * 
 */
enum class StopMode
{
    STEPS, ///< Stop after a fixed number of timesteps
    VOLTAGE ///< Stop when the cell voltage reaches a specified threshold
};


/**
 * @enum GeometryPhase
 * @brief Represents the phase of the geometry in the simulation.
 *
 */
enum class GeometryPhase
{
    SOLID, ///< Solid electrode phase (anode, cathode)
    ELECTROLYTE ///< Electrolyte phase
};

/**
 * @enum BoundarySide
 * @brief Represents the sides of the simulation domain for boundary conditions.
 *
 */
enum class BoundarySide
{
    WEST, ///< West boundary
    EAST, ///< East boundary
    SOUTH, ///< South boundary
    NORTH, ///< North boundary
    BOTTOM, ///< Bottom boundary (for 3D simulations)
    TOP ///< Top boundary (for 3D simulations)
};

/**
 * @enum TIFF_ParticleType
 * @brief Represents the color of the particle in a TIFF image.
 * 
 */
enum class TIFF_ParticleType
{
    BLACK, ///< Particle is black, Electrolye is white 
    WHITE ///< Particle is white, Electrolye is black
};


} // namespace sim
