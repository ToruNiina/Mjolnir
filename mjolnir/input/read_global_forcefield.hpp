#ifndef MJOLNIR_READ_GLOBAL_FORCEFIELD_HPP
#define MJOLNIR_READ_GLOBAL_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_global_interaction.hpp>
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
GlobalForceField<traitsT>
read_global_forcefield(toml::array interactions, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_forcefield(), 0);
    MJOLNIR_LOG_INFO(interactions.size(), " global interactions are found.");

    GlobalForceField<traitsT> gff;
    for(const auto& interaction : interactions)
    {
        if(toml::get<toml::table>(interaction).count("file_name") == 1)
        {
            MJOLNIR_SCOPE(interaction.count("file_name") == 1, 1);

            const auto file_name = toml::find<std::string>(interaction, "file_name");
            if(toml::get<toml::table>(interaction).size() != 1)
            {
                MJOLNIR_LOG_WARN(
                    "[[forcefields.global]] has `file_name` and other keys.");
                MJOLNIR_LOG_WARN(
                    "When `file_name` is provided, other values are ignored "
                    "because those are read from the specified file (",
                    file_name, ").");
            }
            MJOLNIR_LOG_NOTICE("global forcefield is defined in `",
                               input_path, file_name, "`.");

            const auto ff_file = toml::parse(input_path + file_name);
            if(ff_file.count("forcefield") == 1)
            {
                const auto& ff_tab = toml::find(ff_file, "forcefield");
                if(toml::get<toml::table>(ff_tab).count("global") == 1)
                {
                    gff.emplace(read_global_interaction<traitsT>(
                            toml::find(ff_tab, "global")));
                }
            }
            throw_exception<std::runtime_error>("[error] "
                "mjolnir::read_global_forcefield: [forcefield.global] table"
                " should be provided in the file\n --> ", input_path, file_name,
                ".");
        }
        else
        {
            gff.emplace(read_global_interaction<traitsT>(interaction));
        }
    }
    return gff;
}


} // mjolnir
#endif// MJOLNIR_READ_GLOBAL_FORCEFIELD_HPP