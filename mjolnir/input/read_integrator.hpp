#ifndef MJOLNIR_READ_INTEGRATOR
#define MJOLNIR_READ_INTEGRATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletStepper.hpp>
#include <mjolnir/core/UnderdampedLangevinStepper.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
VelocityVerletStepper<traitsT>
read_velocity_verlet_stepper(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_velocity_verlet_stepper(), 0);
    typedef typename traitsT::real_type real_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    return VelocityVerletStepper<traitsT>(delta_t);
}


template<typename traitsT>
UnderdampedLangevinStepper<traitsT>
read_underdamped_langevin_stepper(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_underdamped_langevin_stepper(), 0);
    typedef typename traitsT::real_type real_type;

    const auto seed       = toml::find<std::uint32_t>(simulator, "seed");
    const auto parameters = toml::find<toml::array>(simulator, "parameters");

    std::vector<real_type> gamma(parameters.size());
    for(const auto& params : parameters)
    {
        const auto idx = toml::find<std::size_t>(params, "index");
        const auto  gm = toml::expect<real_type>(params, u8"γ").or_other(
                         toml::expect<real_type>(params, "gamma")).unwrap();
        if(gamma.size() <= idx){gamma.resize(idx+1);}
        gamma.at(idx) = gm;

        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    return UnderdampedLangevinStepper<traitsT>(delta_t, std::move(gamma),
            RandomNumberGenerator<traitsT>(seed));
}

} // mjolnir
#endif// MJOLNIR_READ_INTEGRATOR
