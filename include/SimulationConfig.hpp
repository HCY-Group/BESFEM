/**
 * @file SimulationConfig.hpp
 * @brief Defines the simulation configuration used by BESFEM.
 */

#pragma once

#include "SimTypes.hpp"
#include "Constants.hpp"

#include <string>
#include <vector>

/**
 * @struct SimulationConfig
 * @brief Stores all user-configurable simulation settings.
 *
 * SimulationConfig contains the parameters required to initialize and run a
 * BESFEM simulation, including mesh information, battery configuration,
 * material assignments, initial conditions, and numerical solver settings.
 *
 * The structure is populated from the configuration file and command-line
 * arguments before the simulation begins.
 */
struct SimulationConfig
{
    // -------------------------------------------------------------------------
    // Simulation mode
    // -------------------------------------------------------------------------

    sim::CellMode mode = sim::CellMode::HALF; ///< Cell configuration (half-cell or full-cell).
    sim::Electrode half_electrode = sim::Electrode::ANODE; ///< Active electrode for half-cell simulations.

    // -------------------------------------------------------------------------
    // Input files
    // -------------------------------------------------------------------------

    const char *config_file = "../inputs/run_config.txt"; ///< Simulation configuration file.
    const char *mesh_file   = "../inputs/colored_labels_labels.tif"; ///< Mesh or voxelized geometry file.
    
    // -------------------------------------------------------------------------
    // Discretization
    // -------------------------------------------------------------------------

    int order = Constants::order; ///< Finite element polynomial order.
    int num_timesteps = -1; ///< Number of simulation timesteps.

    bool combine_particle_groups = false; ///< Solve all particle groups as a combined system.

    // -------------------------------------------------------------------------
    // Electrode configuration
    // -------------------------------------------------------------------------

    std::vector<double> init_anode_particles; ///< Initial lithiation of each anode particle.
    std::vector<double> init_cathode_particles; ///< Initial lithiation of each cathode particle.

    std::vector<sim::MaterialType> anode_materials; ///< Material assigned to each anode particle.
    std::vector<sim::MaterialType> cathode_materials; ///< Material assigned to each cathode particle.

    // -------------------------------------------------------------------------
    // Initial conditions
    // -------------------------------------------------------------------------

    double init_CnA = 0.95; ///< Initial anode lithium concentration.
    double init_CnC = 0.30; ///< Initial cathode lithium concentration.
    double init_CnE = -9999; ///< Initial electrolyte concentration.

    double init_BvA = -9999; ///< Initial anode boundary potential (V).
    double init_BvC = -9999; ///< Initial cathode boundary potential (V).
    double init_BvE = -9999; ///< Initial electrolyte potential (V).

    // -------------------------------------------------------------------------
    // Numerical and operating parameters
    // -------------------------------------------------------------------------

    double dh = 5.0e-06; ///< Characteristic mesh spacing (m).
    double gc = 3.3800e-10 * 3.0; ///< Cahn--Hilliard gradient-energy coefficient.
    double dt = 0.001; ///< Simulation timestep.
    double Cr = 1.0; ///< Applied C-rate.
    double Vsr0 = 2.0; ///< Voltage-adjustment rate for constant-current control.

    sim::StopMode stop_mode = sim::StopMode::STEPS; ///< Simulation stopping condition (by steps or voltage).
    double VCut = -1.0; ///< Voltage cutoff for stopping the simulation (V).
    double amr_levels = 0; ///< Number of AMR levels to apply near phase interfaces.

};

/**
 * @brief Parse command-line arguments and configuration file.
 *
 * Reads the simulation configuration file together with any command-line
 * overrides and returns the resulting configuration object.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @return Fully populated simulation configuration.
 */
SimulationConfig ParseSimulationArgs(int argc, char *argv[]);

/**
 * @brief Set up run directly from a config file.
 *
 * Takes the filename of a config file and sets up the run.  
 * Does not accept any additional arguments to override config parameters.
 *
 * @param filename Name of the config file to use.
 * @return Fully populated simulation configuration.
 */
SimulationConfig SimConfigFromFileName(const char *filename);

/**
 * @brief Validate the simulation configuration.
 *
 * Performs consistency checks on the supplied configuration and reports
 * invalid or incompatible options.
 *
 * @param cfg Simulation configuration to validate.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 */
void ValidateConfig(const SimulationConfig &cfg, int argc, char *argv[]);

/**
 * @brief Print the available command-line simulation options.
 */
void PrintAvailableSimulationOptions();
