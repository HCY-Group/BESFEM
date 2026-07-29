#include "../include/SimulationConfig.hpp"
#include "../include/MaterialProperties.hpp"
#include "mfem.hpp"

#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

static std::string Trim(const std::string& s)
{
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";

    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::vector<std::string> SplitString(const std::string& text, char delimiter)
{
    std::vector<std::string> tokens;
    std::stringstream ss(text);
    std::string item;

    while (std::getline(ss, item, delimiter))
    {
        item = Trim(item);
        if (!item.empty())
        {
            tokens.push_back(item);
        }
    }

    return tokens;
}

static std::unordered_map<std::string, std::string>
ReadConfigFile(const std::string& filename)
{
    std::unordered_map<std::string, std::string> data;

    std::ifstream file(filename);
    if (!file)
    {
        mfem::mfem_error(("Could not open config file: " + filename).c_str());
    }

    std::string line;

    while (std::getline(file, line))
    {
        line = Trim(line);

        if (line.empty()) continue;
        if (line[0] == '#') continue;

        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
        {
            line = Trim(line.substr(0, comment_pos));
        }

        const auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key = Trim(line.substr(0, eq_pos));
        std::string val = Trim(line.substr(eq_pos + 1));

        data[key] = val;
    }

    return data;
}

static bool HasKey(
    const std::unordered_map<std::string, std::string>& data,
    const std::string& key)
{
    return data.find(key) != data.end();
}

static std::string GetValue(
    const std::unordered_map<std::string, std::string>& data,
    const std::string& key)
{
    auto it = data.find(key);

    if (it == data.end())
    {
        mfem::mfem_error(("Missing required config key: " + key).c_str());
    }

    return it->second;
}

static bool ParseBool(const std::string& value)
{
    if (value == "true" || value == "1" || value == "yes") return true;
    if (value == "false" || value == "0" || value == "no") return false;

    mfem::mfem_error(("Invalid boolean value: " + value +
                      ". Use true/false.").c_str());

    return false;
}

static sim::MaterialType ParseMaterial(const std::string& name)
{
    if (name == "Graphite" || name == "graphite")
        return sim::MaterialType::Graphite;

    if (name == "NMC" || name == "nmc")
        return sim::MaterialType::NMC;

    if (name == "LFP" || name == "lfp")
        return sim::MaterialType::LFP;

    mfem::mfem_error(("Invalid material: " + name + ". Use Graphite, NMC, or LFP.").c_str());

    return sim::MaterialType::Electrolyte;
}

static std::vector<sim::MaterialType> ParseMaterialList(const std::string& text)
{
    std::vector<sim::MaterialType> materials;

    for (const auto& token : SplitString(text, ','))
    {
        materials.push_back(ParseMaterial(token));
    }

    return materials;
}

static std::vector<double> ParseDoubleList(const std::string& text)
{
    std::vector<double> values;

    for (const auto& token : SplitString(text, ','))
    {
        values.push_back(std::stod(token));
    }

    return values;
}

static sim::StopMode ParseStopMode(const std::string& value)
{
    if (value == "steps")
        return sim::StopMode::STEPS;

    if (value == "voltage")
        return sim::StopMode::VOLTAGE;

    mfem::mfem_error(("Invalid stop_mode: " + value +
                      ". Use steps or voltage.").c_str());

    return sim::StopMode::STEPS;
}

static void ApplyConfigFile(SimulationConfig& cfg)
{
    auto data = ReadConfigFile(cfg.config_file);

    if (HasKey(data, "mode"))
    {
        std::string mode = GetValue(data, "mode");

        if (mode == "half")
            cfg.mode = sim::CellMode::HALF;
        else if (mode == "full")
            cfg.mode = sim::CellMode::FULL;
        else
            mfem::mfem_error("Invalid config mode. Use: half | full.");
    }

    if (HasKey(data, "electrode"))
    {
        std::string electrode = GetValue(data, "electrode");

        if (electrode == "anode")
            cfg.half_electrode = sim::Electrode::ANODE;
        else if (electrode == "cathode")
            cfg.half_electrode = sim::Electrode::CATHODE;
        else
            mfem::mfem_error("Invalid config electrode. Use: anode | cathode.");
    }

    if (HasKey(data, "mesh_file"))
    {
        cfg.mesh_file =
            strdup(GetValue(data, "mesh_file").c_str());
    }

    if (HasKey(data, "num_steps"))
        cfg.num_timesteps = std::stoi(GetValue(data, "num_steps"));

    if (HasKey(data, "combine_particles"))
        cfg.combine_particle_groups = ParseBool(GetValue(data, "combine_particles"));

    if (HasKey(data, "cathode_materials"))
        cfg.cathode_materials = ParseMaterialList(GetValue(data, "cathode_materials"));

    if (HasKey(data, "anode_materials"))
        cfg.anode_materials = ParseMaterialList(GetValue(data, "anode_materials"));

    if (HasKey(data, "init_cathode_particles"))
        cfg.init_cathode_particles = ParseDoubleList(GetValue(data, "init_cathode_particles"));

    if (HasKey(data, "init_anode_particles"))
        cfg.init_anode_particles = ParseDoubleList(GetValue(data, "init_anode_particles"));

    if (HasKey(data, "init_CnE"))
        cfg.init_CnE = std::stod(GetValue(data, "init_CnE"));

    if (HasKey(data, "init_BvA"))
        cfg.init_BvA = std::stod(GetValue(data, "init_BvA"));

    if (HasKey(data, "init_BvC"))
        cfg.init_BvC = std::stod(GetValue(data, "init_BvC"));

    if (HasKey(data, "init_BvE"))
        cfg.init_BvE = std::stod(GetValue(data, "init_BvE"));

    if (HasKey(data, "dh"))
    cfg.dh = std::stod(GetValue(data, "dh"));

    if (HasKey(data, "gc"))
        cfg.gc = std::stod(GetValue(data, "gc"));

    if (HasKey(data, "dt"))
        cfg.dt = std::stod(GetValue(data, "dt"));

    if (HasKey(data, "Cr"))
        cfg.Cr = std::stod(GetValue(data, "Cr"));

    if (HasKey(data, "Vsr0"))
        cfg.Vsr0 = std::stod(GetValue(data, "Vsr0"));

    if (HasKey(data, "stop_mode"))
        cfg.stop_mode = ParseStopMode(GetValue(data, "stop_mode"));

    if (HasKey(data, "VCut"))
        cfg.VCut = std::stod(GetValue(data, "VCut"));

    if (HasKey(data, "amr_levels"))
        cfg.amr_levels = std::stoi(GetValue(data, "amr_levels"));

    if (HasKey(data, "coarsen_factor"))
        cfg.coarsen_factor = std::stoi(GetValue(data, "coarsen_factor"));

    cfg.row_begin =
        std::stoi(GetValue(data, "row_begin"));

    cfg.row_end =
        std::stoi(GetValue(data, "row_end"));

    cfg.column_begin =
        std::stoi(GetValue(data, "column_begin"));

    cfg.column_end =
        std::stoi(GetValue(data, "column_end"));

    if (HasKey(data, "depth_begin"))
        cfg.depth_begin = std::stoi(GetValue(data, "depth_begin"));

    if (HasKey(data, "depth_end"))
        cfg.depth_end = std::stoi(GetValue(data, "depth_end"));

        
}

SimulationConfig ParseSimulationArgs(int argc, char *argv[])
{
    SimulationConfig cfg;

    const char* config_file = "../inputs/run_config.txt";
    bool list_options = false;

    for (int i = 1; i < argc; ++i)
    {
        if ((std::strcmp(argv[i], "-cfg") == 0 ||
            std::strcmp(argv[i], "--config") == 0) &&
            i + 1 < argc)
        {
            config_file = argv[i + 1];
            ++i;
        }
        else if (std::strcmp(argv[i], "-opts") == 0 ||
                std::strcmp(argv[i], "--list-options") == 0)
        {
            list_options = true;
        }
    }

    cfg.config_file = config_file;

    if (list_options)
    {
        PrintAvailableSimulationOptions();
        std::exit(0);
    }

    ApplyConfigFile(cfg);

    const char* mode =
        (cfg.mode == sim::CellMode::FULL) ? "full" : "half";

    const char* half_elec =
        (cfg.half_electrode == sim::Electrode::CATHODE) ? "cathode" : "anode";

    const char* stop_mode =
        (cfg.stop_mode == sim::StopMode::STEPS) ? "steps" : "voltage";

    mfem::OptionsParser args(argc, argv);

    args.AddOption(&cfg.config_file,
                   "-cfg", "--config",
                   "User-editable BESFEM config file.");

    args.AddOption(&list_options,
                   "-opts", "--list-options",
                   "-no-opts", "--no-list-options",
                   "Print available BESFEM option choices.");

    args.AddOption(&cfg.mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&cfg.order, "-o", "--order", "Finite element polynomial degree.");
    args.AddOption(&mode, "-mode", "--mode", "Cell mode: half | full.");
    args.AddOption(&half_elec, "-elec", "--electrode", "HALF mode only: anode | cathode.");
    args.AddOption(&cfg.combine_particle_groups, "-combine", "--combine-particles", "-separate", "--separate-particles", "Combine all particle groups into one.");
    args.AddOption(&stop_mode, "-stop", "--stop-mode", "Stopping mode: steps | voltage.");

    args.ParseCheck();

    if (std::strcmp(mode, "half") == 0)
        cfg.mode = sim::CellMode::HALF;
    else if (std::strcmp(mode, "full") == 0)
        cfg.mode = sim::CellMode::FULL;
    else
        mfem::mfem_error("Invalid -mode. Use: half | full.");

    if (std::strcmp(half_elec, "anode") == 0)
        cfg.half_electrode = sim::Electrode::ANODE;
    else if (std::strcmp(half_elec, "cathode") == 0)
        cfg.half_electrode = sim::Electrode::CATHODE;
    else
        mfem::mfem_error("Invalid -elec. Use: anode | cathode.");

    if (std::strcmp(stop_mode, "steps") == 0)
        cfg.stop_mode = sim::StopMode::STEPS;
    else if (std::strcmp(stop_mode, "voltage") == 0)
        cfg.stop_mode = sim::StopMode::VOLTAGE;
    else
        mfem::mfem_error("Invalid stop mode. Use: steps | voltage.");

    return cfg;
}

static bool IsCathodeMaterial(sim::MaterialType m)
{
    return m == sim::MaterialType::NMC ||
           m == sim::MaterialType::LFP;
}

static bool IsAnodeMaterial(sim::MaterialType m)
{
    return m == sim::MaterialType::Graphite;
}
static void CheckCathodeInitialBoundaryFromOCV(const SimulationConfig& cfg)
{
    const double tolerance = 0.3;

    for (size_t i = 0; i < cfg.cathode_materials.size(); i++)
    {
        sim::MaterialType material = cfg.cathode_materials[i];
        double x = cfg.init_cathode_particles[i];

        double expected_ocv = MaterialProperties::OCV(material, x);

        double init_cathode_ocv = cfg.init_BvC - cfg.init_BvE;

        if (std::abs(init_cathode_ocv - expected_ocv) > tolerance)
        {
            std::stringstream ss;
            ss << "Initial cathode potential difference is inconsistent with "
                  "the OCV curve for cathode particle "
               << i << ". For x = " << x
               << ", expected init_BvC - init_BvE to be about "
               << expected_ocv
               << ", but got "
               << init_cathode_ocv
               << ". Suggested init_BvC is about "
               << cfg.init_BvE + expected_ocv << ".";

            mfem::mfem_error(ss.str().c_str());
        }
    }
}
static void CheckAnodeInitialBoundaryFromOCV(const SimulationConfig& cfg)
{
    const double tolerance = 0.1;

    for (size_t i = 0; i < cfg.anode_materials.size(); i++)
    {
        sim::MaterialType material = cfg.anode_materials[i];
        double x = cfg.init_anode_particles[i];

        double expected_ocv = MaterialProperties::OCV(material, x);

        double initial_anode_ocv = cfg.init_BvA - cfg.init_BvE;

        bool match_pos = std::abs(initial_anode_ocv - expected_ocv) <= tolerance;
        bool match_neg = std::abs(initial_anode_ocv + expected_ocv) <= tolerance;

        if (!match_pos && !match_neg)
        {
            double suggested_pos_BvA = cfg.init_BvE + expected_ocv;
            double suggested_neg_BvA = cfg.init_BvE - expected_ocv;
            std::stringstream ss;
            ss << "Initial anode potential difference is inconsistent with "
              "the OCV curve for anode particle "
            << i << ".\n"
            << "For x = " << x
            << ", expected init_BvA - init_BvE to be about +"
            << expected_ocv << " or -" << expected_ocv
            << ", but got " << initial_anode_ocv << ".\n"
            << "Suggested init_BvA values are about "
            << suggested_pos_BvA
            << " or "
            << suggested_neg_BvA
            << ".";
        }
    }
}
static void CheckParticleStoichiometry(
    const std::vector<double>& values,
    const std::string& name)
{
    for (size_t i = 0; i < values.size(); i++)
    {
        double x = values[i];

        if (x < 0.0 || x > 1.0)
        {
            std::stringstream ss;
            ss << name << "[" << i << "] = " << x
               << " is invalid. Initial particle stoichiometry "
               << "must be between 0 and 1.";

            mfem::mfem_error(ss.str().c_str());
        }
    }
}

void ValidateConfig(const SimulationConfig &cfg, int argc, char *argv[])
{
    if (cfg.mode != sim::CellMode::HALF)
    {
        mfem::mfem_error(
            "Only HALF-CATHODE mode is currently implemented.");
    }
    
    const bool cathode = cfg.half_electrode == sim::Electrode::CATHODE;

    

    if (!cfg.mesh_file)
        mfem::mfem_error("mesh_file cannot be empty.");

    if (cfg.stop_mode == sim::StopMode::STEPS)
    {
        if (cfg.num_timesteps <= 0)
        {
            mfem::mfem_error(
                "stop_mode=steps requires num_steps > 0.");
        }
    }
    else if (cfg.stop_mode == sim::StopMode::VOLTAGE)
    {
        if (cfg.VCut <= 0.0)
        {
            mfem::mfem_error(
                "stop_mode=voltage requires VCut > 0.");
        }
    }

    if (cfg.amr_levels < 0)
    {
        mfem::mfem_error(
            "amr_levels must be greater than or equal to 0.");
    }

    if (cfg.coarsen_factor < 1)
    {
        mfem::mfem_error(
            "coarsen_factor must be an integer greater than or equal to 1.");
    }

    // Coarsening is only permitted when AMR is enabled.
    if (cfg.amr_levels == 0 && cfg.coarsen_factor != 1)
    {
        mfem::mfem_error(
            "coarsen_factor must equal 1 when amr_levels = 0. "
            "Initial mesh coarsening is only allowed when AMR is enabled.");
    }

    // ------------------------------------------------------------
    // TIFF cropping
    //
    // For each direction:
    //   -1, -1  = use the entire TIFF extent
    //   >=0     = use the specified [begin, end) crop
    // ------------------------------------------------------------

    const bool use_all_rows =
        cfg.row_begin == -1 &&
        cfg.row_end == -1;

    const bool crop_rows =
        cfg.row_begin >= 0 &&
        cfg.row_end >= 0;

    const bool use_all_columns =
        cfg.column_begin == -1 &&
        cfg.column_end == -1;

    const bool crop_columns =
        cfg.column_begin >= 0 &&
        cfg.column_end >= 0;

    // Rows must either both be -1 or both be nonnegative.
    if (!use_all_rows && !crop_rows)
    {
        mfem::mfem_error(
            "row_begin and row_end must either both be -1 "
            "to use all rows, or both be nonnegative crop bounds.");
    }

    // Columns must either both be -1 or both be nonnegative.
    if (!use_all_columns && !crop_columns)
    {
        mfem::mfem_error(
            "column_begin and column_end must either both be -1 "
            "to use all columns, or both be nonnegative crop bounds.");
    }

    // Check ordering only when an actual crop was requested.
    if (crop_rows && cfg.row_end <= cfg.row_begin)
    {
        mfem::mfem_error(
            "Invalid row crop: row_end must be greater than row_begin.");
    }

    if (crop_columns && cfg.column_end <= cfg.column_begin)
    {
        mfem::mfem_error(
            "Invalid column crop: column_end must be greater than column_begin.");
    }

    // Require at least two selected points in a cropped direction.
    if (crop_rows && cfg.row_end - cfg.row_begin < 2)
    {
        mfem::mfem_error(
            "The row crop must contain at least two TIFF points.");
    }

    if (crop_columns && cfg.column_end - cfg.column_begin < 2)
    {
        mfem::mfem_error(
            "The column crop must contain at least two TIFF points.");
    }

    if (cfg.init_BvE == -9999.0)
    {
        mfem::mfem_error(
            "Missing init_BvE in run_config.txt.");
    }
    if (cfg.init_CnE == -9999.0)
    {
        mfem::mfem_error(
            "Missing init_CnE in run_config.txt.");
    }
    if (cfg.dt <= 0.0)
    mfem::mfem_error("dt must be positive.");
    
    if (cfg.dh <= 0.0)
        mfem::mfem_error("dh must be positive.");
    
    if (cfg.gc <= 0.0)
        mfem::mfem_error("gc must be positive.");
    
    if (cfg.Cr <= 0.0)
        mfem::mfem_error("Cr must be positive.");
    if (cfg.stop_mode == sim::StopMode::STEPS)
    {
        if (cfg.num_timesteps <= 0)
        {
            mfem::mfem_error(
                "stop_mode=steps requires num_steps > 0.");
        }
    }
    else if (cfg.stop_mode == sim::StopMode::VOLTAGE)
    {
        if (cfg.VCut <= 0.0)
        {
            mfem::mfem_error(
                "stop_mode=voltage requires VCut > 0.");
        }
    }

    if (cfg.mode == sim::CellMode::FULL)
    {
        if (cfg.cathode_materials.empty())
            mfem::mfem_error("cathode_materials cannot be empty, please list these in your config file.");

        if (cfg.anode_materials.empty())
            mfem::mfem_error("anode_materials cannot be empty, please list these in your config file.");

        if (cfg.init_cathode_particles.empty())
            mfem::mfem_error("init_cathode_particles cannot be empty, please list these in your config file.");

        if (cfg.init_anode_particles.empty())
            mfem::mfem_error("init_anode_particles cannot be empty, please list these in your config file.");
    }
    else 
    {
        if (cathode)
        {
    
            if (cfg.cathode_materials.empty())
                mfem::mfem_error("cathode_materials cannot be empty.");
    
            if (cfg.init_cathode_particles.empty())
                mfem::mfem_error("init_cathode_particles cannot be empty.");

            if (cfg.cathode_materials.size() != cfg.init_cathode_particles.size())
            {
                mfem::mfem_error(
                    "cathode_materials and init_cathode_particles must have the same length.");
            }
             if (cfg.init_BvC == -9999.0)
            {
                mfem::mfem_error(
                    "Missing init_BvC. HALF-CATHODE mode requires init_BvC in run_config.txt.");
            }
            
            for (auto m : cfg.cathode_materials)
            {
                if (!IsCathodeMaterial(m))
                {
                    mfem::mfem_error(
                        "Invalid cathode_materials: cathode can only use NMC or LFP. "
                        "Graphite is an anode material.");
                }
            }
            CheckCathodeInitialBoundaryFromOCV(cfg);
             CheckParticleStoichiometry(
                    cfg.init_cathode_particles,
                    "init_cathode_particles");
            
        }
        else
        {
            if (cfg.anode_materials.empty())
                mfem::mfem_error("anode_materials cannot be empty.");
    
            if (cfg.init_anode_particles.empty())
                mfem::mfem_error("init_anode_particles cannot be empty.");
            if (cfg.anode_materials.size() != cfg.init_anode_particles.size())
            {
                mfem::mfem_error(
                    "anode_materials and init_anode_particles must have the same length.");
            }
            // if (cfg.anode_materials.size() != 3)
            // {
            //     std::stringstream ss;
            //     ss << "Invalid anode particle group count. "
            //        << "BESFEM currently expects exactly 3 anode psi/material groups, but you provided "
            //        << cfg.anode_materials.size()
            //        << ". Example: anode_materials = Graphite,Graphite,Graphite ";     
            //     mfem::mfem_error(ss.str().c_str());
            // }
            // if (cfg.init_anode_particles.size() != 3)
            // {
            //     std::stringstream ss;
            //     ss << "Invalid anode particle group count. "
            //        << "BESFEM currently expects exactly 3 anode particles groups, but you provided "
            //        << cfg.init_anode_particles.size()
            //        << ". Example: init_anode_particles = 0.02,0.02,0.02.";
            
            //     mfem::mfem_error(ss.str().c_str());
            // }
             if (cfg.init_BvA == -9999.0)
            {
                mfem::mfem_error(
                    "Missing init_BvA. HALF-CATHODE mode requires init_BvA in run_config.txt.");
            }
            for (auto m : cfg.anode_materials)
            {
                if (!IsAnodeMaterial(m))
                {
                    mfem::mfem_error(
                        "Invalid anode_materials: anode can only use graphite.");
                }
            }
            CheckAnodeInitialBoundaryFromOCV(cfg);
            CheckParticleStoichiometry(
                    cfg.init_anode_particles,
                    "init_anode_particles");
    
        }
    }

    
}

void PrintAvailableSimulationOptions()
{
    std::cout << "\nAvailable BESFEM command-line choices to edit in ../inputs/run_config.txt:\n\n";

    std::cout << "  Run config:\n";
    std::cout << "    -cfg ../inputs/run_config.txt\n\n";

    std::cout << "  Cell modes:\n";
    std::cout << "    half\n\n";

    std::cout << "  Half-cell electrodes:\n";
    std::cout << "    anode\n";
    std::cout << "    cathode\n\n";

    std::cout << "  Materials:\n";
    std::cout << "    cathode_materials: NMC, LFP\n";
    std::cout << "    anode_materials:   Graphite\n\n";

    std::cout << "  Mixed cathode example:\n";
    std::cout << "    cathode_materials = LFP,LFP,NMC\n\n";

    std::cout << "  Mesh types:\n";
    std::cout << "    ml   MATLAB mesh\n";
    std::cout << "    v    voxel mesh\n\n";

    std::cout << "  Particle grouping:\n";
    std::cout << "    combine_particles = true\n";
    std::cout << "    combine_particles = false\n\n";

    std::cout << "  Initial particle concentrations:\n";
    std::cout << "    init_cathode_particles = 0.2,0.15,0.10\n";
    std::cout << "    init_anode_particles = 0.8,0.85,0.90\n\n";

    std::cout << "  Number of timesteps:\n";
    std::cout << "    num_steps = 1000\n\n";

    std::cout << " Initial boundary values:\n";
    std::cout << "    init_BvA = -0.01\n";
    std::cout << "    init_BvC = 3.30\n";
    std::cout << "    init_BvE = -0.10\n\n";

    std::cout << " Simulation Parameters:\n";
    std::cout << "    timestep size dt = 0.001\n";
    std::cout << "    element size dh = 2.0e-6\n";
    std::cout << "    Charge Rate Cr = 1.0\n";
    std::cout << "    Voltage Scanning Rate Vsr0 = 0.009\n";
    std::cout << "    Cahn Hilliard Gradient Coef gc = 1.014e-9\n\n";

    std::cout << "  Stopping criteria:\n";
    std::cout << "    stop_mode = steps\n";
    std::cout << "    stop_mode = voltage\n";
    std::cout << "    VCut = 3.2\n\n";

    std::cout << "  Example command:\n";
    std::cout << "    mpirun -np 4 ./battery_simulation -cfg ../inputs/run_config.txt\n\n";
}