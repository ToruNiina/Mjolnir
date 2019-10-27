#include <mjolnir/omp_nondet/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<OpenMPNonDeterministicSimulatorTraits<float,  UnlimitedBoundary>>>;
} // mjolnir
