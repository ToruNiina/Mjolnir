#include <mjolnir/forcefield/external/PositionRestraintInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PositionRestraintInteraction<SimulatorTraits<double, UnlimitedBoundary       >, HarmonicPotential<double>>;
template class PositionRestraintInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, HarmonicPotential<float> >;
template class PositionRestraintInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
template class PositionRestraintInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;
} // mjolnir
