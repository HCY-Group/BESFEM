#pragma once
#include "SimTypes.hpp"
#include "Constants.hpp"
#include <string>
#include <vector>

/**
 * @file SimulationConfig.hpp
 * @brief Defines the simulation configuration structure and argument parsing utilities.
 *
 * Provides a lightweight configuration container for selecting the cell mode,
 * electrode type, mesh files, phase-field files, finite element order, and
 * timestep count. Also declares helper routines for parsing and validating
 * command-line arguments.
 */

/**
 * @struct SimulationConfig
 * @brief Holds all user-specified input parameters for a BESFEM run.
 *
 * This small POD-style structure stores simulation parameters such as:
 * - cell mode (half-cell or full-cell),
 * - active electrode in half-cell mode,
 * - mesh / distance-function file paths,
 * - mesh type indicator used by Domain_Parameters,
 * - finite element polynomial order,
 * - total number of time steps.
 *
 * Default values are taken from constants in `Constants.hpp`.
 */
struct SimulationConfig {
    sim::CellMode mode = sim::CellMode::HALF; ///< Simulation mode (HALF or FULL).
    sim::Electrode half_electrode = sim::Electrode::ANODE; ///< Electrode used in HALF mode.
    
    const char* mesh_file   = Constants::mesh_file;   ///< Path to mesh file.
    const char* dsF_file_A  = Constants::dsF_file_A;  ///< Distance function file for anode.
    const char* dsF_file_C  = Constants::dsF_file_C;  ///< Distance function file for cathode.

    const char* mesh_type = nullptr; ///< Mesh type specifier ("ml" or "v").
    
    int order = Constants::order; ///< Finite element polynomial order.
    int num_timesteps = 1000;     ///< Number of global time steps for the simulation.


    std::vector<double> init_anode_particles;
    std::vector<double> init_cathode_particles;
};

/**
 * @brief Parse command-line arguments into a SimulationConfig structure.
 *
 * Expected supported arguments include:
 * - `--full` or `--half`
 * - `--anode`, `--cathode`, or `--both`
 * - `--mesh <file>`
 * - `--dsA <file>` and `--dsC <file>`
 * - `--order <p>`
 * - `--steps <N>`
 *
 * Missing arguments are replaced by defaults from `Constants.hpp`.
 *
 * @param argc Standard argument count.
 * @param argv Standard argument vector.
 * @return A fully populated SimulationConfig structure.
 */
SimulationConfig ParseSimulationArgs(int argc, char *argv[]);

/**
 * @brief Validate that the chosen configuration is logically consistent.
 *
 * Checks include:
 * - HALF mode must specify either ANODE or CATHODE.
 * - FULL mode cannot use `half_electrode == BOTH` explicitly.
 * - Required distance-function files must be present when using FULL mode.
 * - Mesh type must be non-null.
 *
 * If validation fails, this function prints an error message and terminates.
 *
 * @param cfg  Simulation configuration to validate.
 * @param argc Argument count (for error messages).
 * @param argv Argument vector (for error messages).
 */
void ValidateConfig(const SimulationConfig &cfg, int argc, char *argv[]);

