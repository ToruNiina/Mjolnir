#ifndef MJOLNIR_READ_OBSERVER
#define MJOLNIR_READ_OBSERVER
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
Observer<traitsT>
read_observer(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_observer(), 0);

    const auto& general =
        get_toml_value<toml::Table>(data, "general", "<root>");

    std::string path =
        get_toml_value<std::string>(general, "output_path", "[general]");
    if(path.back() != '/') {path += '/';} //XXX assuming posix

    const std::string prefix = get_toml_value<std::string>(
            general, "output_prefix", "[general]");
    MJOLNIR_LOG_INFO("path   = ", path);
    MJOLNIR_LOG_INFO("prefix = ", prefix);

    return Observer<traitsT>(path + prefix);
}



}
#endif// MJOLNIR_READ_OBSERVER
