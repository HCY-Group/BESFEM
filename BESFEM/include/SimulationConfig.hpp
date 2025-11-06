#pragma once
#include "SimTypes.hpp"
#include "../inputs/Constants.hpp"
#include <string>

struct SimulationConfig {
    sim::CellMode mode = sim::CellMode::HALF;
    sim::Electrode half_electrode = sim::Electrode::ANODE;
    const char* mesh_file = Constants::mesh_file;
    const char* dsF_file_A = Constants::dsF_file_A;
    const char* dsF_file_C = Constants::dsF_file_C;
    const char* mesh_type = nullptr;
    int order = Constants::order;
    int num_timesteps = 1000;
};

// Parse command-line args and return a filled config
SimulationConfig ParseSimulationArgs(int argc, char *argv[]);

// Validate that configuration is consistent with mode and electrode
void ValidateConfig(const SimulationConfig &cfg, int argc, char *argv[]);
