#include <mjolnir/core/GJFNPTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GJFNPTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GJFNPTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
