#include "../include/MaterialProperties.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#ifndef BESFEM_MATERIALS_DIR
#define BESFEM_MATERIALS_DIR "inputs/materials"
#endif

namespace MaterialProperties
{
    namespace
    {
        const std::vector<std::string> material_names = {
            "Graphite", "LFP", "NMC", "Carbon", "Silicon", "Electrolyte"};
        const std::vector<sim::MaterialType> material_types = {
            sim::MaterialType::Graphite, sim::MaterialType::LFP,
            sim::MaterialType::NMC, sim::MaterialType::Carbon,
            sim::MaterialType::Silicon, sim::MaterialType::Electrolyte};
        const std::vector<std::string> property_names = {
            "ocv", "chemical_potential", "exchange_current_density",
            "diffusivity", "mobility", "conductivity", "site_density"};

        // The same index in each vector belongs to the same property.
        std::vector<sim::MaterialType> materials;
        std::vector<std::string> names;
        std::vector<double> constants;
        std::vector<std::vector<double>> concentrations;
        std::vector<std::vector<double>> values;
        bool initialized = false;

        void ReadFile(const std::filesystem::path& path, size_t index)
        {
            std::ifstream file(path);
            std::string line;
            while (std::getline(file, line))
            {
                // Ignore comments and blank lines.
                line = line.substr(0, line.find('#'));
                std::istringstream row(line);
                double first;
                double second;
                if (row >> first)
                {
                    if (row >> second)
                    {
                        concentrations[index].push_back(first);
                        values[index].push_back(second);
                    }
                    else
                    {
                        constants[index] = first;
                    }
                }
            }
        }

        double Evaluate(sim::MaterialType material, const std::string& name, double concentration)
        {
            if (!initialized)
            {
                Configure({}, "");
            }

            for (size_t i = 0; i < names.size(); i++)
            {
                if (materials[i] != material || names[i] != name)
                {
                    continue;
                }
                if (concentrations[i].empty())
                {
                    return constants[i];
                }

                // Use the endpoint value outside the table.
                if (concentration <= concentrations[i].front())
                {
                    return values[i].front();
                }
                if (concentration >= concentrations[i].back())
                {
                    return values[i].back();
                }

                // Find the two surrounding points and interpolate between them.
                for (size_t j = 1; j < concentrations[i].size(); j++)
                {
                    if (concentration <= concentrations[i][j])
                    {
                        double x1 = concentrations[i][j - 1];
                        double x2 = concentrations[i][j];
                        double y1 = values[i][j - 1];
                        double y2 = values[i][j];
                        double fraction = (concentration - x1) / (x2 - x1);
                        return y1 + fraction * (y2 - y1);
                    }
                }
            }
            return 0.0;
        }
    }

    void Configure(const std::unordered_map<std::string, std::string>& settings,
                   const std::string& config_file)
    {
        std::filesystem::path base = std::filesystem::current_path();
        if (!config_file.empty())
        {
            base = std::filesystem::absolute(config_file).parent_path();
        }
        std::filesystem::path root = BESFEM_MATERIALS_DIR;
        if (settings.count("materials_dir") > 0)
        {
            root = base / settings.at("materials_dir");
        }

        materials.clear();
        names.clear();
        constants.clear();
        concentrations.clear();
        values.clear();

        for (size_t i = 0; i < material_names.size(); i++)
        {
            for (size_t j = 0; j < property_names.size(); j++)
            {
                size_t index = names.size();
                materials.push_back(material_types[i]);
                names.push_back(property_names[j]);
                constants.push_back(0.0);
                concentrations.push_back({});
                values.push_back({});

                std::string key = "material." + material_names[i] + "." + property_names[j];
                if (settings.count(key) > 0)
                {
                    std::istringstream setting(settings.at(key));
                    std::string kind;
                    setting >> kind;
                    if (kind == "constant")
                    {
                        setting >> constants[index];
                    }
                    else if (kind == "table")
                    {
                        std::string filename;
                        setting >> std::quoted(filename);
                        ReadFile(base / filename, index);
                    }
                }
                else
                {
                    std::filesystem::path path = root / material_names[i] / (property_names[j] + ".txt");
                    ReadFile(path, index);
                }
            }
        }
        initialized = true;
    }

    double OCV(sim::MaterialType material, double concentration)
    {
        return Evaluate(material, "ocv", concentration);
    }

    double ChemicalPotential(sim::MaterialType material, double concentration)
    {
        return Evaluate(material, "chemical_potential", concentration);
    }

    double ExchangeCurrentDensity(sim::MaterialType material, double concentration)
    {
        return Evaluate(material, "exchange_current_density", concentration);
    }

    double Diffusivity(sim::MaterialType material, double concentration)
    {
        return Evaluate(material, "diffusivity", concentration);
    }

    double Mobility(sim::MaterialType material, double concentration)
    {
        return Evaluate(material, "mobility", concentration);
    }

    double Conductivity(sim::MaterialType material, double concentration)
    {
        return Evaluate(material, "conductivity", concentration);
    }

    double SiteDensity(sim::MaterialType material)
    {
        return Evaluate(material, "site_density", 0.0);
    }
}
