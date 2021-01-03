#include <mjolnir/omp/GJFNPTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GJFNPTLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GJFNPTLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
