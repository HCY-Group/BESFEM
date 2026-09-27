#include "../include/MaterialProperties.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifndef BESFEM_MATERIALS_DIR
#define BESFEM_MATERIALS_DIR "inputs/materials"
#endif

namespace MaterialProperties
{
    // Helpers and stored data below are private to this source file.
    namespace
    {
        using Material = sim::MaterialType;
        struct MaterialInfo
        {
            std::string name;
            Material type;
        };

        const std::vector<MaterialInfo> materials = {
            {"Graphite", Material::Graphite}, {"LFP", Material::LFP},
            {"NMC", Material::NMC}, {"Carbon", Material::Carbon},
            {"Silicon", Material::Silicon}, {"Electrolyte", Material::Electrolyte}};
        const std::vector<std::string> property_names = {
            "ocv", "chemical_potential", "exchange_current_density",
            "diffusivity", "mobility", "conductivity", "site_density"};

        // Only list properties implemented by the original models. Missing
        // physics is not silently replaced by a zero-valued property.
        bool HasDefault(Material material, const std::string& name)
        {
            if (material == Material::Electrolyte)
            {
                return name == "diffusivity";
            }
            if (name == "mobility")
            {
                return material == Material::Graphite || material == Material::LFP;
            }
            if (name == "diffusivity")
            {
                return material != Material::Graphite;
            }
            return true;
        }

        // Each property holds either one constant or a table of values.
        // An empty concentration vector means this is a constant property.
        struct Property1D
        {
            Material material;
            std::string name;
            double constant = 0.0;
            std::vector<double> concentrations;
            std::vector<double> values;
            bool error_outside = false;

            double Evaluate(double concentration) const
            {
                if (!std::isfinite(concentration))
                {
                    throw std::runtime_error("nonfinite concentration");
                }
                if (concentrations.empty())
                {
                    return constant;
                }

                // With the "error" option, reject values outside the table.
                // Otherwise, clamp them to the nearest endpoint below.
                if (error_outside && (concentration < concentrations.front() || concentration > concentrations.back()))
                {
                    throw std::runtime_error("concentration outside table range [" +
                        std::to_string(concentrations.front()) + ", " +
                        std::to_string(concentrations.back()) + "]");
                }
                if (concentration <= concentrations.front())
                {
                    return values.front();
                }
                if (concentration >= concentrations.back())
                {
                    return values.back();
                }

                // Find the first table concentration at or above the requested one.
                // lower_bound uses binary search, which is fast for large tables.
                auto upper_point = std::lower_bound(concentrations.begin(), concentrations.end(), concentration);
                size_t upper_index = upper_point - concentrations.begin();
                size_t lower_index = upper_index - 1;

                // Linear interpolation: move the same fraction between the values
                // as the requested concentration lies between the concentrations.
                double fraction = (concentration - concentrations[lower_index]) / (concentrations[upper_index] - concentrations[lower_index]);
                double value_change = values[upper_index] - values[lower_index];
                return values[lower_index] + fraction * value_change;
            }
        };

        // One entry per loaded property (for example, NMC's OCV).
        std::vector<Property1D> database;
        bool initialized = false;

        void ValidateValue(const std::string& name, double value)
        {
            if (!std::isfinite(value))
            {
                throw std::runtime_error("value must be finite");
            }
            if (name == "site_density" && value <= 0)
            {
                throw std::runtime_error("site density must be positive");
            }
            bool requires_nonnegative_value =
                name == "diffusivity" || name == "mobility" ||
                name == "conductivity" || name == "exchange_current_density";
            if (requires_nonnegative_value && value < 0)
            {
                throw std::runtime_error("value must be nonnegative");
            }
        }

        Property1D ReadFile(Property1D property, const std::filesystem::path& path)
        {
            std::ifstream input(path);
            if (!input)
            {
                throw std::runtime_error("cannot open " + path.string());
            }
            std::string line;
            size_t line_number = 0;
            bool scalar = false;
            while (std::getline(input, line))
            {
                ++line_number;
                // Remove comments, then skip blank or whitespace-only lines.
                size_t comment_start = line.find('#');
                line = line.substr(0, comment_start);
                std::istringstream row(line);
                row >> std::ws;
                if (row.eof())
                {
                    continue;
                }
                try
                {
                    double concentration;
                    double value;
                    if (scalar)
                    {
                        throw std::runtime_error("scalar file must contain exactly one number");
                    }
                    if (!(row >> concentration))
                    {
                        throw std::runtime_error("expected a number");
                    }
                    row >> std::ws;
                    // A first row with just one number defines a constant.
                    // Otherwise every data row must contain concentration/value.
                    if (row.eof() && property.concentrations.empty())
                    {
                        ValidateValue(property.name, concentration);
                        property.constant = concentration;
                        scalar = true;
                        continue;
                    }
                    std::string extra;
                    if (!(row >> value) || (row >> extra))
                    {
                        throw std::runtime_error("expected two numbers");
                    }
                    if (property.name == "site_density")
                    {
                        throw std::runtime_error("site density file must contain one scalar");
                    }
                    if (!std::isfinite(concentration) || (!property.concentrations.empty() && concentration <= property.concentrations.back()))
                    {
                        throw std::runtime_error("concentrations must be finite and strictly increasing");
                    }
                    if (concentration < 0 || (property.material != Material::Electrolyte && concentration > 1))
                    {
                        throw std::runtime_error("concentration outside physical domain");
                    }
                    ValidateValue(property.name, value);
                    property.concentrations.push_back(concentration);
                    property.values.push_back(value);
                }
                catch (const std::runtime_error& error)
                {
                    throw std::runtime_error(path.string() + ":" + std::to_string(line_number) + ": " + error.what());
                }
            }
            if (input.bad())
            {
                throw std::runtime_error("error reading " + path.string());
            }
            if (!scalar && property.concentrations.size() < 2)
            {
                throw std::runtime_error(path.string() + ": need a scalar or at least two concentration/value rows");
            }
            return property;
        }

        Property1D ReadOverride(Material material, const std::string& name,
                                const std::string& setting, const std::filesystem::path& base)
        {
            Property1D property;
            property.material = material;
            property.name = name;

            std::istringstream input(setting);
            std::string kind, extra;
            input >> kind;
            if (kind == "constant")
            {
                if (!(input >> property.constant) || (input >> extra))
                {
                    throw std::runtime_error("expected constant <number>");
                }
                ValidateValue(name, property.constant);
                return property;
            }
            if (kind != "table")
            {
                throw std::runtime_error("expected constant or table");
            }

            std::string filename;
            std::string policy = "clamp";
            if (!(input >> std::quoted(filename)))
            {
                throw std::runtime_error("expected table <path> [clamp|error]");
            }
            input >> policy;
            if ((policy != "clamp" && policy != "error") || (input >> extra))
            {
                throw std::runtime_error("expected table <path> [clamp|error]");
            }
            property.error_outside = policy == "error";
            return ReadFile(property, base / filename);
        }

        double Evaluate(Material material, const std::string& name, double concentration)
        {
            if (!initialized)
            {
                Configure({}, "");
            }
            for (const Property1D& property : database)
            {
                if (property.material == material && property.name == name)
                {
                    return property.Evaluate(concentration);
                }
            }
            throw std::runtime_error("No property defined for " + name);
        }
    }

    void Configure(const std::unordered_map<std::string, std::string>& settings,
                   const std::string& config_file)
    {
        std::string config_path = config_file;
        if (config_path.empty())
        {
            config_path = ".";
        }
        std::filesystem::path base = std::filesystem::absolute(config_path).parent_path();
        std::filesystem::path root = BESFEM_MATERIALS_DIR;
        if (settings.count("materials_dir") > 0)
        {
            root = base / settings.at("materials_dir");
        }

        // Keep the old data until all new files have loaded successfully.
        std::vector<Property1D> loaded;
        size_t overrides_loaded = 0;
        for (const MaterialInfo& material : materials)
        {
            for (const std::string& name : property_names)
            {
                std::string setting_name = "material." + material.name + "." + name;
                bool overridden = settings.count(setting_name) > 0;
                if (material.type == Material::Electrolyte && name != "diffusivity")
                {
                    continue;
                }
                try
                {
                    if (overridden)
                    {
                        loaded.push_back(ReadOverride(material.type, name, settings.at(setting_name), base));
                        ++overrides_loaded;
                    }
                    else if (HasDefault(material.type, name))
                    {
                        Property1D property;
                        property.material = material.type;
                        property.name = name;
                        loaded.push_back(ReadFile(property, root / material.name / (name + ".txt")));
                    }
                }
                catch (const std::runtime_error& error)
                {
                    throw std::runtime_error(setting_name + ": " + error.what());
                }
            }
        }

        // Catch misspelled or unsupported material settings instead of ignoring them.
        size_t overrides_supplied = 0;
        for (const auto& setting : settings)
        {
            if (setting.first.compare(0, 9, "material.") == 0)
            {
                ++overrides_supplied;
            }
        }
        if (overrides_loaded != overrides_supplied)
        {
            throw std::runtime_error("Unknown or unsupported material property setting");
        }
        database.swap(loaded);
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
