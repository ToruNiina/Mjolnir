#include <mjolnir/forcefield/external/ExternalDistanceInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;


template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;


template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

} // mjolnir
