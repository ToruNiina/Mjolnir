#include <mjolnir/omp_nondet/BAOABLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class BAOABLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class BAOABLangevinIntegrator<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
