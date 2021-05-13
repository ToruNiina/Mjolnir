#include <mjolnir/forcefield/global/GlobalPairWCAInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        WCAPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        WCAPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, WCAPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, WCAPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

} // mjolnir
