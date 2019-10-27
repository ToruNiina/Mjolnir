#include <mjolnir/omp_nondet/RandomNumberGenerator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class RandomNumberGenerator<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>>;
template class RandomNumberGenerator<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>>;
template class RandomNumberGenerator<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class RandomNumberGenerator<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
