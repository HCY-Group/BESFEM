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

    // Half-cell geometry.
    const char *mesh_file = nullptr;

    // Full-cell geometries.
    const char *anode_mesh_file = nullptr;
    const char *cathode_mesh_file = nullptr;
    
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
    double gc = -1.0; ///< Cahn--Hilliard gradient-energy coefficient.
    double dt = 0.001; ///< Simulation timestep.
    double Cr = 1.0; ///< Applied C-rate.
    double Vsr0 = 2.0; ///< Voltage-adjustment rate for constant-current control.

    sim::StopMode stop_mode = sim::StopMode::STEPS; ///< Simulation stopping condition (by steps or voltage).
    double VCut = -1.0; ///< Voltage cutoff for stopping the simulation (V).
    double amr_levels = 0; ///< Number of AMR levels to apply near phase interfaces.

    // -------------------------------------------------------------------------
    // AMR parameters & geometry cropping
    // -------------------------------------------------------------------------

    // Initial structured-grid coarsening factor.
    // 1 means no coarsening, 2 keeps every second TIFF node, etc.
    int coarsen_factor = 1; ///< Initial structured-grid coarsening factor.

    int row_begin = -1; ///< Beginning row for geometry cropping.
    int row_end = -1;   ///< Ending row for geometry cropping.
    int column_begin = -1; ///< Beginning column for geometry cropping.
    int column_end = -1;   ///< Ending column for geometry cropping.

    // Optional for future 3D use.
    int depth_begin = 0; ///< Beginning depth for geometry cropping.
    int depth_end = 1;   ///< Ending depth for geometry cropping.

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