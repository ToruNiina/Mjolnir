#ifndef MJOLNIR_TEST_MAKE_EMPTY_INPUT_HPP
#define MJOLNIR_TEST_MAKE_EMPTY_INPUT_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/string.hpp>

namespace mjolnir
{
namespace test
{

// make a minimal input file that contains no particle, no forcefield.
// It contains only parameters that are needed to pass read_input() functions.
inline toml::value make_empty_input()
{
    using namespace mjolnir::literals::string_literals;

    toml::table output;
    output["prefix"] = "nothing"_s;
    output["format"] = "xyz"_s;
    output["path"]   = "./"_s;

    toml::table files;
    files["output"]  = std::move(output);

    toml::table units;
    units["length"] = "angstrom"_s;
    units["energy"] = "kcal/mol"_s;

    // Mjolnir does not consider running simulation without simulator
    // (that does not make sense!) . Here, temporary set MD simulator with
    // Newtonian dynamics.
    toml::table integrator;
    integrator["type"] = "VelocityVerlet"_s;
    toml::table simulator;
    simulator["type"]          = "MolecularDynamics"_s;
    simulator["integrator"]    = std::move(integrator);
    simulator["precision"]     = "double"_s;
    simulator["boundary_type"] = "Unlimited"_s;
    simulator["total_step"]    = 1;
    simulator["save_step"]     = 1;
    simulator["delta_t"]       = 0.1;

    // Mjolnir requires a system to simulate, even if it has no particle.
    toml::table system;
    system["attributes"]     = toml::table{{"temperature"_s, toml::value(300.0)}};
    system["boundary_shape"] = toml::table{}; // no boundary
    system["particles"]      = toml::array{}; // no particle
    toml::array systems{std::move(system)};

    // "empty forcefields" completely make sense because essentially any kind of
    // forcefields are a sum of several potential terms. Forcefield that has
    // zero term means an ideal system having no interaction.
    toml::array forcefields{toml::table{/* nothing! */}};

    toml::table input;
    input["files"]       = std::move(files);
    input["units"]       = std::move(units);
    input["simulator"]   = std::move(simulator);
    input["systems"]     = std::move(systems);
    input["forcefields"] = std::move(forcefields);

    return toml::value(input);
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_MAKE_EMPTY_INPUT_HPP
