#include <mjolnir/cuda/ObserverContainer.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ObserverContainer<CUDASimulatorTraits<double, UnlimitedBoundary>       >;
template class ObserverContainer<CUDASimulatorTraits<float,  UnlimitedBoundary>       >;
template class ObserverContainer<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ObserverContainer<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
