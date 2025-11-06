#include "../include/SimulationConfig.hpp"
#include "mfem.hpp"
#include <cstring>
#include <filesystem>
#include <fstream>

SimulationConfig ParseSimulationArgs(int argc, char *argv[])
{
    SimulationConfig cfg;
    const char* mode = "half";
    const char* half_elec = "anode";

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&cfg.mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&cfg.dsF_file_C, "-dC", "--cathode-distance", "Cathode distance file.");
    args.AddOption(&cfg.dsF_file_A, "-dA", "--anode-distance", "Anode distance file.");
    args.AddOption(&cfg.order, "-o", "--order", "Finite element polynomial degree.");
    args.AddOption(&cfg.mesh_type, "-t", "--type", "Mesh type (r, c, v, d).");
    args.AddOption(&cfg.num_timesteps, "-n", "--num-steps", "Number of timesteps.");
    args.AddOption(&mode, "-mode", "--mode", "Cell mode: half | full.");
    args.AddOption(&half_elec, "-elec", "--electrode", "HALF mode only: anode | cathode.");
    args.ParseCheck();

    cfg.mode = (std::strcmp(mode, "full") == 0)
               ? sim::CellMode::FULL : sim::CellMode::HALF;
    cfg.half_electrode = (std::strcmp(half_elec, "cathode") == 0)
               ? sim::Electrode::CATHODE : sim::Electrode::ANODE;

    return cfg;
}

void ValidateConfig(const SimulationConfig &cfg, int argc, char *argv[])
{
    auto flag_present = [&](std::initializer_list<const char*> names)
    {
        for (int i = 1; i < argc; ++i)
            for (auto n : names)
                if (std::strcmp(argv[i], n) == 0) return true;
        return false;
    };

    const bool used_dA = flag_present({"-dA", "--anode-distance"});
    const bool used_dC = flag_present({"-dC", "--cathode-distance"});

    if (cfg.mode == sim::CellMode::FULL) {
        if (!used_dA || !used_dC)
            mfem::mfem_error("FULL mode requires both -dA and -dC.");
    } else {
        bool cathode = (cfg.half_electrode == sim::Electrode::CATHODE);
        if (cathode && !used_dC)
            mfem::mfem_error("HALF-CATHODE requires -dC <file>.");
        if (!cathode && !used_dA)
            mfem::mfem_error("HALF-ANODE requires -dA <file>.");
    }
}
