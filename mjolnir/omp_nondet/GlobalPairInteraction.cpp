#include <mjolnir/omp_nondet/GlobalPairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// D-H
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
