#include "cece/physics/cece_speciation_config.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "conf/config.hpp"

namespace cece {

SpeciationConfig SpeciationConfigLoader::Load(const std::string& mechanism_path, const std::string& mapping_path, const std::string& dataset) const {
    if (!std::filesystem::exists(mechanism_path)) throw std::runtime_error("Mechanism file not found: " + mechanism_path);
    if (!std::filesystem::exists(mapping_path)) throw std::runtime_error("Mapping file not found: " + mapping_path);

    conf::Config mechanism = conf::Config::from_file(mechanism_path);
    conf::Config mapping = conf::Config::from_file(mapping_path);
    SpeciationConfig config = ParseMechanism(mechanism.root());
    ParseMapping(mapping.root(), config, dataset);
    Validate(config);
    return config;
}

SpeciationConfig SpeciationConfigLoader::ParseMechanism(const conf::Value& node) const {
    SpeciationConfig config;
    if (!node["name"]) throw std::invalid_argument("Mechanism file missing required 'name' key");
    config.mechanism_name = node["name"].as_string();

    conf::Value species = node["species"];
    if (!species || species.kind() != conf::Node_Kind::Sequence) throw std::invalid_argument("Mechanism file missing required 'species' list");
    for (std::size_t i = 0; i < species.size(); ++i) {
        conf::Value entry = species[i];
        if (!entry["name"]) throw std::invalid_argument("Mechanism species entry " + std::to_string(i) + " missing required 'name' field");
        if (!entry["molecular weight [kg mol-1]"]) {
            throw std::invalid_argument("Mechanism species entry " + std::to_string(i) + " missing required 'molecular weight [kg mol-1]' field");
        }

        MechanismSpecies sp;
        sp.name = entry["name"].as<std::string>();

        // Convert from kg/mol (MICM convention) to g/mol
        double mw_kg_per_mol = entry["molecular weight [kg mol-1]"].as<double>();
        if (mw_kg_per_mol <= 0.0) {
            throw std::invalid_argument("Mechanism species '" + sp.name + "' has non-positive molecular weight: " + std::to_string(mw_kg_per_mol));
        }
        sp.molecular_weight = mw_kg_per_mol * 1000.0;  // kg/mol -> g/mol

        // Aerosol species declare themselves explicitly and carry size-bin properties.
        if (entry["is_aerosol"] && entry["is_aerosol"].as<bool>()) {
            sp.is_aerosol = true;
            for (const char* key : {"density [kg m-3]", "lower_radius [um]", "upper_radius [um]"}) {
                if (!entry[key]) {
                    throw std::invalid_argument("Aerosol species '" + sp.name + "' missing required '" + key + "' field");
                }
            }
            sp.density = entry["density [kg m-3]"].as<double>();
            sp.lower_radius = entry["lower_radius [um]"].as<double>();
            sp.upper_radius = entry["upper_radius [um]"].as<double>();
            sp.effective_radius =
                entry["effective_radius [um]"] ? entry["effective_radius [um]"].as<double>() : 0.5 * (sp.lower_radius + sp.upper_radius);

            if (sp.density <= 0.0) {
                throw std::invalid_argument("Aerosol species '" + sp.name + "' has non-positive density");
            }
            if (sp.lower_radius <= 0.0 || sp.upper_radius <= 0.0) {
                throw std::invalid_argument("Aerosol species '" + sp.name + "' has non-positive radius bounds");
            }
            if (sp.lower_radius >= sp.upper_radius) {
                throw std::invalid_argument("Aerosol species '" + sp.name + "' requires lower_radius < upper_radius");
            }
        }

        config.species.push_back(sp);
    }
    return config;
}

void SpeciationConfigLoader::ParseMapping(const conf::Value& node, SpeciationConfig& config, const std::string& dataset) const {
    if (!node["mechanism"]) throw std::invalid_argument("Mapping file missing required 'mechanism' key");
    conf::Value datasets = node["datasets"];
    if (!datasets || datasets.kind() != conf::Node_Kind::Map) throw std::invalid_argument("Mapping file missing required 'datasets' section");
    conf::Value selected = datasets[dataset];
    if (!selected) throw std::invalid_argument("Requested dataset '" + dataset + "' not found in mapping file");
    if (selected.kind() != conf::Node_Kind::Map) throw std::invalid_argument("Dataset '" + dataset + "' is not a map");
    config.dataset_name = dataset;

    for (const auto& mechanism_species : selected.keys()) {
        conf::Value class_map = selected[mechanism_species];
        if (class_map.kind() != conf::Node_Kind::Map) {
            throw std::invalid_argument("Mechanism species '" + mechanism_species + "' in dataset '" + dataset +
                                        "' is not a map of emission classes");
        }

        // Iterate emission class → scale factor pairs
        for (auto class_it = class_map.begin(); class_it != class_map.end(); ++class_it) {
            const auto& key_node = class_it->first;

            // Skip null nodes
            if (key_node.Type() == YAML::NodeType::Null || !key_node.IsDefined()) {
                continue;
            }

            std::string class_name = key_node.Scalar();
            if (class_name.empty()) {
                try {
                    class_name = key_node.as<std::string>();
                } catch (...) {
                    class_name = "";
                }
            }

            // Handle yaml-cpp YAML 1.1 boolean interpretation of "NO"
            if (class_name == "false" || class_name == "no") {
                class_name = "NO";
            }

            EmissionClass ec;
            bool resolved = StringToEmissionClass(class_name, ec);
            if (!resolved) {
                std::string upper_name = class_name;
                std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(), ::toupper);
                resolved = StringToEmissionClass(upper_name, ec);
            }
            if (!resolved) {
                // Aerosol identity source: the class names an aerosol species (size
                // bin) defined in the mechanism rather than a gas emission class.
                bool is_aerosol_source = false;
                for (const auto& sp : config.species) {
                    if (sp.name == class_name && sp.is_aerosol) {
                        is_aerosol_source = true;
                        break;
                    }
                }
                if (!is_aerosol_source) {
                    throw std::invalid_argument("Invalid emission class '" + class_name + "' for mechanism species '" + mechanism_species +
                                                "' in dataset '" + dataset + "'");
                }
                ec = EmissionClass::COUNT;  // sentinel: aerosol identity source, not a gas class
            }
            double scale_factor = class_map[class_name].as_double();
            if (scale_factor <= 0.0) throw std::invalid_argument("Non-positive scale factor for emission class '" + class_name + "'");
            config.mappings.push_back({mechanism_species, emission_class, scale_factor});
        }
    }
}

void SpeciationConfigLoader::Validate(const SpeciationConfig& config) const {
    std::unordered_set<std::string> species_names;
    for (const auto& species : config.species) species_names.insert(species.name);
    std::vector<std::string> unknown;
    for (const auto& mapping : config.mappings)
        if (!species_names.contains(mapping.mechanism_species)) unknown.push_back(mapping.mechanism_species);
    if (!unknown.empty()) {
        std::sort(unknown.begin(), unknown.end());
        unknown.erase(std::unique(unknown.begin(), unknown.end()), unknown.end());
        std::ostringstream message;
        message << "Mapping references unknown mechanism species not in mechanism file: ";
        for (std::size_t i = 0; i < unknown.size(); ++i) message << (i ? ", " : "") << "'" << unknown[i] << "'";
        throw std::invalid_argument(message.str());
    }
    for (const auto& mapping : config.mappings) {
        if (mechanism_species_names.find(mapping.mechanism_species) == mechanism_species_names.end()) {
            unknown_species.push_back(mapping.mechanism_species);
        }
    }

    if (!unknown_species.empty()) {
        std::sort(unknown_species.begin(), unknown_species.end());
        unknown_species.erase(std::unique(unknown_species.begin(), unknown_species.end()), unknown_species.end());
        std::ostringstream oss;
        oss << "Mapping references unknown mechanism species not in mechanism file: ";
        for (std::size_t i = 0; i < unknown_species.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << "'" << unknown_species[i] << "'";
        }
        throw std::invalid_argument(oss.str());
    }

    // Emission class validation is already done during ParseMapping (StringToEmissionClass),
    // but we double-check here for configs built programmatically
    for (const auto& mapping : config.mappings) {
        // COUNT is the aerosol-identity sentinel (no gas emission class).
        if (mapping.emission_class == EmissionClass::COUNT) continue;
        int ec_idx = static_cast<int>(mapping.emission_class);
        if (ec_idx < 0 || ec_idx >= static_cast<int>(EmissionClass::COUNT)) {
            throw std::invalid_argument("Invalid emission class index " + std::to_string(ec_idx) + " in mapping for mechanism species '" +
                                        mapping.mechanism_species + "'");
        }
    }
}

std::string SpeciationConfigLoader::ToYaml(const SpeciationConfig& config) {
    // Direct stream formatting avoids reintroducing yaml-cpp emitter dependency
    std::ostringstream out;
    out << "name: \"" << config.mechanism_name << "\"\nspecies:\n";
    for (const auto& species : config.species) {
        out << "  - name: \"" << species.name << "\"\n";
        out << "    molecular weight [kg mol-1]: " << species.molecular_weight / 1000.0 << "\n";
    }
    out << "mechanism: \"" << config.mechanism_name << "\"\ndatasets:\n  \"" << config.dataset_name << "\":\n";
    std::unordered_map<std::string, std::vector<const SpeciationMapping*>> grouped;
    for (const auto& mapping : config.mappings) grouped[mapping.mechanism_species].push_back(&mapping);
    for (const auto& [species, mappings] : grouped) {
        out << "    \"" << species << "\":\n";
        for (const auto* mapping : mappings)
            out << "      \"" << EmissionClassToString(mapping->emission_class) << "\": " << mapping->scale_factor << "\n";
    }
    return out.str();
}

}  // namespace cece
