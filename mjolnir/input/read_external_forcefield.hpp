#ifndef MJOLNIR_READ_EXTERNAL_FORCEFIELD_HPP
#define MJOLNIR_READ_EXTERNAL_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_external_interaction.hpp>
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
ExternalForceField<traitsT>
read_external_forcefield(toml::array interactions, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_external_forcefield(), 0);
    MJOLNIR_LOG_INFO(interactions.size(), " external interactions are found.");

    ExternalForceField<traitsT> eff;
    for(const auto& interaction : interactions)
    {
        if(toml::get<toml::table>(interaction).count("file_name") == 1)
        {
            MJOLNIR_SCOPE(interaction.count("file_name") == 1, 1);

            const auto file_name = toml::find<std::string>(interaction, "file_name");
            if(toml::get<toml::table>(interaction).size() != 1)
            {
                MJOLNIR_LOG_WARN(
                    "[[forcefields.external]] has `file_name` and other keys.");
                MJOLNIR_LOG_WARN(
                    "When `file_name` is provided, other values are ignored "
                    "because those are read from the specified file (",
                    file_name, ").");
            }
            MJOLNIR_LOG_NOTICE("external forcefield is defined in `",
                               input_path, file_name, "`.");

            const auto ff_file = toml::parse(input_path + file_name);
            if(ff_file.count("forcefield") == 1)
            {
                const auto& ff_tab = toml::find(ff_file, "forcefield");
                if(toml::get<toml::table>(ff_tab).count("external") == 1)
                {
                    eff.emplace(read_external_interaction<traitsT>(
                            toml::find(ff_tab, "external")));
                }
            }
            throw_exception<std::runtime_error>("[error] "
                "mjolnir::read_external_forcefield: [forcefield.external] table"
                " should be provided in the file\n --> ", input_path, file_name,
                ".");
        }
        else
        {
            eff.emplace(read_external_interaction<traitsT>(interaction));
        }
    }
    return eff;
}

} // mjolnir
#endif// MJOLNIR_READ_EXTERNAL_FORCEFIELD_HPP