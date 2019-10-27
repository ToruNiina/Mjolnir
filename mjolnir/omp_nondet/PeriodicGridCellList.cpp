#include <mjolnir/omp_nondet/PeriodicGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
