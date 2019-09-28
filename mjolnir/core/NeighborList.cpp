#include <mjolnir/core/NeighborList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class NeighborList<empty_t, std::uint32_t>;
template class NeighborList<float,  std::uint32_t>;
template class NeighborList<double, std::uint32_t>;
template class NeighborList<std::pair<float , float >, std::uint32_t>;
template class NeighborList<std::pair<double, double>, std::uint32_t>;
} // mjolnir
