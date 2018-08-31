#ifndef MJOLNIR_NEIGHBOR_LIST_HPP
#define MJOLNIR_NEIGHBOR_LIST_HPP
#include <mjolnir/util/range.hpp>
#include <mjolnir/util/empty.hpp>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

//! a container for verlet-list, cell-list, and other spatial indexing methods.
//! the objectives are
//! - store indices of possible partners that has interaction.
//! - store pre-calculated parameters. e.g. epsilon_ij = sqrt(eps_i * eps_j)
//!   for Lennard-Jones takes relatively large computational cost to obtain.
//!   calculate it for all the possible pairs to make force calculation fast.

namespace mjolnir
{

// google with "empty base optimization(EBO)".
// by using EBO, it compress the object size when `paramT` is a empty class.
// without this, empty parameter type consumes 8 byte (on x86_64) because
// std::size_t requires 8-byte alignment.

namespace detail
{
// In order to handle `double` or other non-class object, neighbor_element need
// another layer that checks `std::is_empty`. If empty, it inherits it and EBO
// removes the overhead. If not empty(like `double`), it keeps storage for the
// value. `neighbor_elmenet_impl` stores value if `std::true_type` is passed.
template<typename paramT, bool IsEmpty>
struct neighbor_element_impl;

template<typename paramT>
struct neighbor_element_impl<paramT, true> : private paramT // for EBO
{
    using parameter_type = paramT;

    neighbor_element_impl(const parameter_type& p) {}
    neighbor_element_impl(parameter_type&& p)      {}

    parameter_type&       parameter()       noexcept {return *this;}
    parameter_type const& parameter() const noexcept {return *this;}
};

template<typename paramT>
struct neighbor_element_impl<paramT, false>
{
    using parameter_type = paramT;

    neighbor_element_impl(const parameter_type& p) : param(p) {}
    neighbor_element_impl(parameter_type&& p)      : param(std::move(p)) {}

    parameter_type&       parameter()       noexcept {return param;}
    parameter_type const& parameter() const noexcept {return param;}

    paramT param;
};

} // detail

template<typename paramT>
struct neighbor_element
: private detail::neighbor_element_impl<paramT, std::is_empty<paramT>::value>
{
    using base_type =
        detail::neighbor_element_impl<paramT, std::is_empty<paramT>::value>;
    using parameter_type = typename base_type::parameter_type;

    neighbor_element(std::size_t idx, std::size_t jdx, const paramT& p)
        : base_type(p), i(idx), j(jdx)
    {}
    neighbor_element(std::size_t idx, std::size_t jdx, paramT&& p)
        : base_type(std::move(p)), i(idx), j(jdx)
    {}

    neighbor_element()  = default;
    ~neighbor_element() = default;
    neighbor_element(const neighbor_element&) = default;
    neighbor_element(neighbor_element&&)      = default;
    neighbor_element& operator=(const neighbor_element&) = default;
    neighbor_element& operator=(neighbor_element&&)      = default;

    parameter_type&       parameter()       noexcept {return base_type::parameter();}
    parameter_type const& parameter() const noexcept {return base_type::parameter();}

    std::size_t i, j;
};

template<typename paramT>
inline bool operator==(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return (lhs.i == rhs.i) && (lhs.j == rhs.j);
}
template<typename paramT>
inline bool operator!=(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return !(lhs == rhs);
}
template<typename paramT>
inline bool operator<(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return (lhs.i == rhs.i) ? (lhs.j < rhs.j) : (lhs.i < rhs.i);
}
template<typename paramT>
inline bool operator<=(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return (lhs == rhs) || (lhs < rhs);
}
template<typename paramT>
inline bool operator>(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return !(lhs <= rhs);
}
template<typename paramT>
inline bool operator>=(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return !(lhs < rhs);
}

template<typename parameterT>
class NeighborList
{
  public:
    using parameter_type = parameterT;
    using neighbor_type  = neighbor_element<parameter_type>;
    using container_type = std::vector<neighbor_type>;
    using iterator       = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

  public:

    NeighborList() = default;
    ~NeighborList() = default;
    NeighborList(const NeighborList&) = default;
    NeighborList(NeighborList&&)      = default;
    NeighborList& operator=(const NeighborList&) = default;
    NeighborList& operator=(NeighborList&&)      = default;

    void clear() {this->neighbors_.clear();}

    // if number of particles inside of cutoff length could be known, it might
    // increase runtime speed a bit
    void reserve(std::size_t Nparticle, std::size_t Nneighbor)
    {
        this->neighbors_.reserve(Nparticle * Nneighbor);
    }

    void push_back(const neighbor_type& n)
    {
        this->neighbors_.push_back(n);
    }
    void push_back(neighbor_type&& n)
    {
        this->neighbors_.push_back(std::move(n));
    }

    template<typename ... Ts>
    void emplace_back(Ts&& ... args)
    {
        this->neighbors_.emplace_back(std::forward<Ts>(args)...);
    }

    template<typename Iterator>
    void append(Iterator first, Iterator last)
    {
        static_assert(std::is_same<neighbor_type,
            typename std::iterator_traits<Iterator>::value_type>::value,
            "iterator value_type must be the same as neighbor_type");

        std::copy(first, last, std::back_inserter(this->neighbors_));
        return ;
    }

    const_iterator  begin() const noexcept {return neighbors_.cbegin();}
    const_iterator  end()   const noexcept {return neighbors_.cend();}
    const_iterator cbegin() const noexcept {return neighbors_.cbegin();}
    const_iterator cend()   const noexcept {return neighbors_.cend();}

  private:
    container_type neighbors_;
};

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
