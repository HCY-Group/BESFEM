#include "../include/MaterialProperties.hpp"
#include "../include/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifndef BESFEM_MATERIALS_DIR
#define BESFEM_MATERIALS_DIR "inputs/materials"
#endif

namespace MaterialProperties
{
    namespace
    {
        using Material = sim::MaterialType;
        using Key = std::pair<Material, std::string>;
        const std::map<std::string, Material> materials = {
            {"Graphite", Material::Graphite}, {"LFP", Material::LFP},
            {"NMC", Material::NMC}, {"Carbon", Material::Carbon},
            {"Silicon", Material::Silicon}, {"Electrolyte", Material::Electrolyte}};
        const std::vector<std::string> property_names = {
            "ocv", "chemical_potential", "exchange_current_density",
            "diffusivity", "mobility", "conductivity", "site_density"};

        // Only list properties implemented by the original models. Missing
        // physics is not silently replaced by a zero-valued property.
        bool HasDefault(Material m, const std::string& name)
        {
            if (m == Material::Electrolyte) return name == "diffusivity";
            if (name == "mobility") return m == Material::Graphite || m == Material::LFP;
            if (name == "diffusivity") return m != Material::Graphite;
            return true;
        }

        std::string Label(const Key& key)
        {
            for (const auto& material : materials)
                if (material.second == key.first) return material.first + "." + key.second;
            return "unknown material." + key.second;
        }

        struct Property1D
        {
            double constant = 0.0;
            std::vector<double> x, y;
            bool error_outside = false;
            std::string source;

            double Evaluate(double c) const
            {
                if (!std::isfinite(c)) throw std::runtime_error("nonfinite concentration");
                if (x.empty()) return constant;
                if (error_outside && (c < x.front() || c > x.back()))
                    throw std::runtime_error("concentration outside table range [" +
                        std::to_string(x.front()) + ", " + std::to_string(x.back()) + "]");
                if (c <= x.front()) return y.front();
                if (c >= x.back()) return y.back();
                const auto i = std::lower_bound(x.begin(), x.end(), c) - x.begin();
                return y[i-1] + (c-x[i-1])/(x[i]-x[i-1])*(y[i]-y[i-1]);
            }
        };

        std::map<Key, Property1D> database;
        bool initialized = false;

        void ValidateValue(const Key& key, double value)
        {
            const auto& name = key.second;
            if (!std::isfinite(value)) throw std::runtime_error("value must be finite");
            if (name == "site_density" && value <= 0)
                throw std::runtime_error("site density must be positive");
            if ((name == "diffusivity" || name == "mobility" ||
                 name == "conductivity" || name == "exchange_current_density") && value < 0)
                throw std::runtime_error("value must be nonnegative");
        }

        Property1D ReadProperty(const Key& key, const std::string& specification,
                                const std::filesystem::path& base)
        {
            Property1D property;
            property.source = specification;
            std::istringstream spec(specification);
            std::string kind, extra;
            spec >> kind;
            if (kind == "constant")
            {
                if (!(spec >> property.constant) || (spec >> extra))
                    throw std::runtime_error("expected constant <number>");
                ValidateValue(key, property.constant);
                return property;
            }
            if (kind != "table") throw std::runtime_error("expected constant or table");
            std::string filename, policy;
            if (!(spec >> std::quoted(filename))) throw std::runtime_error("expected table <path> [clamp|error]");
            if (spec >> policy)
            {
                if (policy != "clamp" && policy != "error")
                    throw std::runtime_error("expected clamp or error");
                property.error_outside = policy == "error";
            }
            if (spec >> extra) throw std::runtime_error("unexpected table options");
            const auto path = base / filename;
            property.source = path.string();
            std::ifstream input(path);
            if (!input) throw std::runtime_error("cannot open " + path.string());
            std::string line;
            size_t line_number = 0;
            bool scalar = false;
            while (std::getline(input, line))
            {
                ++line_number;
                line = line.substr(0, line.find('#'));
                std::istringstream row(line);
                row >> std::ws;
                if (row.eof()) continue;
                try
                {
                    double x, y;
                    if (scalar) throw std::runtime_error("scalar file must contain exactly one number");
                    if (!(row >> x)) throw std::runtime_error("expected a number");
                    row >> std::ws;
                    if (row.eof() && property.x.empty())
                    {
                        ValidateValue(key, x);
                        property.constant = x;
                        scalar = true;
                        continue;
                    }
                    if (!(row >> y) || (row >> extra)) throw std::runtime_error("expected two numbers");
                    if (key.second == "site_density")
                        throw std::runtime_error("site density file must contain one scalar");
                    if (!std::isfinite(x) || (!property.x.empty() && x <= property.x.back()))
                        throw std::runtime_error("concentrations must be finite and strictly increasing");
                    if (x < 0 || (key.first != Material::Electrolyte && x > 1))
                        throw std::runtime_error("concentration outside physical domain");
                    ValidateValue(key, y);
                    property.x.push_back(x);
                    property.y.push_back(y);
                }
                catch (const std::runtime_error& error)
                {
                    throw std::runtime_error(path.string() + ":" +
                        std::to_string(line_number) + ": " + error.what());
                }
            }
            if (input.bad()) throw std::runtime_error("error reading " + path.string());
            if (!scalar && property.x.size() < 2)
                throw std::runtime_error(path.string() + ": need a scalar or at least two concentration/value rows");
            return property;
        }

        double Evaluate(Material material, const std::string& name, double c)
        {
            if (!initialized) Configure({}, "");
            const Key key{material, name};
            const auto it = database.find(key);
            if (it == database.end()) throw std::runtime_error("No property defined for " + Label(key));
            try { return it->second.Evaluate(c); }
            catch (const std::runtime_error& error)
            {
                throw std::runtime_error(Label(key) + " (" + it->second.source +
                    ") at concentration " + std::to_string(c) + ": " + error.what());
            }
        }
    }

    void Configure(const std::unordered_map<std::string, std::string>& values,
                   const std::string& config_file)
    {
        const auto base = std::filesystem::absolute(config_file.empty() ? "." : config_file).parent_path();
        auto root = std::filesystem::path(BESFEM_MATERIALS_DIR);
        const auto directory = values.find("materials_dir");
        if (directory != values.end()) root = base / directory->second;
        root = std::filesystem::absolute(root);

        // Assemble the final sources before reading: an override replaces its
        // default file entirely, even when that default file is absent.
        std::map<Key, std::string> specifications;
        for (const auto& material : materials)
            for (const auto& name : property_names)
                if (HasDefault(material.second, name))
                {
                    std::ostringstream spec;
                    spec << "table " << std::quoted((root / material.first / (name + ".txt")).string());
                    specifications[{material.second, name}] = spec.str();
                }

        std::map<Key, std::string> explicit_overrides;
        for (const auto& entry : values)
        {
            if (entry.first.compare(0, 9, "material.") != 0) continue;
            const auto dot = entry.first.find('.', 9);
            const auto material = materials.find(entry.first.substr(9, dot - 9));
            const auto name = dot == std::string::npos ? "" : entry.first.substr(dot + 1);
            if (material == materials.end()) throw std::runtime_error(entry.first + ": unknown material");
            if (std::find(property_names.begin(), property_names.end(), name) == property_names.end() ||
                (name == "chp_value" && material->second != Material::LFP))
                throw std::runtime_error(entry.first + ": unknown property");
            if (material->second == Material::Electrolyte && name != "diffusivity")
                throw std::runtime_error(entry.first + ": only electrolyte diffusivity is supported");
            const Key key{material->second, name};
            explicit_overrides[key] = entry.second;
            specifications[key] = entry.second;
        }

        // Preserve the previously supported OCV override dependency. Default
        // chemical potentials themselves are now independent, editable files.
        for (const auto& material : materials)
            if (material.second != Material::Graphite &&
                explicit_overrides.count({material.second, "ocv"}) &&
                !explicit_overrides.count({material.second, "chemical_potential"}))
                specifications.erase({material.second, "chemical_potential"});

        std::map<Key, Property1D> next;
        for (const auto& specification : specifications)
        {
            try { next[specification.first] = ReadProperty(specification.first, specification.second, base); }
            catch (const std::runtime_error& error)
            {
                throw std::runtime_error(Label(specification.first) + ": " + error.what());
            }
        }
        for (const auto& material : materials)
        {
            const auto m = material.second;
            if (m == Material::Graphite || !explicit_overrides.count({m, "ocv"}) ||
                explicit_overrides.count({m, "chemical_potential"})) continue;
            auto mu = next.at({m, "ocv"});
            const double scale = m == Material::Carbon ? -1.0 : -Constants::Frd;
            mu.constant *= scale;
            ValidateValue({m, "chemical_potential"}, mu.constant);
            for (auto& value : mu.y)
            {
                value *= scale;
                ValidateValue({m, "chemical_potential"}, value);
            }
            mu.source = "derived from " + mu.source;
            next[{m, "chemical_potential"}] = std::move(mu);
        }
        database.swap(next);
        initialized = true;
    }

    double OCV(sim::MaterialType m, double c) { return Evaluate(m, "ocv", c); }
    double ChemicalPotential(sim::MaterialType m, double c) { return Evaluate(m, "chemical_potential", c); }
    double ExchangeCurrentDensity(sim::MaterialType m, double c) { return Evaluate(m, "exchange_current_density", c); }
    double Diffusivity(sim::MaterialType m, double c) { return Evaluate(m, "diffusivity", c); }
    double Mobility(sim::MaterialType m, double c) { return Evaluate(m, "mobility", c); }
    double Conductivity(sim::MaterialType m, double c) { return Evaluate(m, "conductivity", c); }
    double SiteDensity(sim::MaterialType m) { return Evaluate(m, "site_density", 0.0); }
    double LFP_ChpValue(double c) { return Evaluate(sim::MaterialType::LFP, "chp_value", c); }
}
