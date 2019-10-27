#include <mjolnir/omp_nondet/GlobalPairLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
