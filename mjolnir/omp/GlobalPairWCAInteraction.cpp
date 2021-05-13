#include <mjolnir/omp/GlobalPairWCAInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       , WCAPotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       , WCAPotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, WCAPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, WCAPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
