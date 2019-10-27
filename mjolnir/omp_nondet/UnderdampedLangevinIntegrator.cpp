#include <mjolnir/omp_nondet/UnderdampedLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnderdampedLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>>;
template class UnderdampedLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>>;
template class UnderdampedLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class UnderdampedLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
