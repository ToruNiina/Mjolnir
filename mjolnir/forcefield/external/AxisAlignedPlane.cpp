#include <mjolnir/forcefield/external/AxisAlignedPlane.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class AxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary>       , XAxis<SimulatorTraits<double, UnlimitedBoundary>       , PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary>       , XAxis<SimulatorTraits<float,  UnlimitedBoundary>       , PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>, XAxis<SimulatorTraits<double, CuboidalPeriodicBoundary>, PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>, XAxis<SimulatorTraits<float,  CuboidalPeriodicBoundary>, PositiveDirection>>;

template class AxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary>       , XAxis<SimulatorTraits<double, UnlimitedBoundary>       , NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary>       , XAxis<SimulatorTraits<float,  UnlimitedBoundary>       , NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>, XAxis<SimulatorTraits<double, CuboidalPeriodicBoundary>, NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>, XAxis<SimulatorTraits<float,  CuboidalPeriodicBoundary>, NegativeDirection>>;

template class AxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary>       , YAxis<SimulatorTraits<double, UnlimitedBoundary>       , PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary>       , YAxis<SimulatorTraits<float,  UnlimitedBoundary>       , PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>, YAxis<SimulatorTraits<double, CuboidalPeriodicBoundary>, PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>, YAxis<SimulatorTraits<float,  CuboidalPeriodicBoundary>, PositiveDirection>>;

template class AxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary>       , YAxis<SimulatorTraits<double, UnlimitedBoundary>       , NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary>       , YAxis<SimulatorTraits<float,  UnlimitedBoundary>       , NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>, YAxis<SimulatorTraits<double, CuboidalPeriodicBoundary>, NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>, YAxis<SimulatorTraits<float,  CuboidalPeriodicBoundary>, NegativeDirection>>;

template class AxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary>       , ZAxis<SimulatorTraits<double, UnlimitedBoundary>       , PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary>       , ZAxis<SimulatorTraits<float,  UnlimitedBoundary>       , PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>, ZAxis<SimulatorTraits<double, CuboidalPeriodicBoundary>, PositiveDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ZAxis<SimulatorTraits<float,  CuboidalPeriodicBoundary>, PositiveDirection>>;

template class AxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary>       , ZAxis<SimulatorTraits<double, UnlimitedBoundary>       , NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary>       , ZAxis<SimulatorTraits<float,  UnlimitedBoundary>       , NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>, ZAxis<SimulatorTraits<double, CuboidalPeriodicBoundary>, NegativeDirection>>;
template class AxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ZAxis<SimulatorTraits<float,  CuboidalPeriodicBoundary>, NegativeDirection>>;

} // mjolnir
