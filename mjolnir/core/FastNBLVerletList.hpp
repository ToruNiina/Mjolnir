#ifndef MJOLNIR_CORE_UNLIMITED_FASTNBL_VERLET_LIST_HPP
#define MJOLNIR_CORE_UNLIMITED_FASTNBL_VERLET_LIST_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <algorithm>
#include <limits>

namespace mjolnir
{

// Verlet list based FastNBL implementation. This class is for UnlimitedBoundary case.
// FastNBL is introduced in the following paper.
// - Kun Li, Shigang Li, Shan Huang, Yifeng Chen, Yunquan Zhang (2019) J. Supercomput.
template<typename traitsT, typename parameterT>
class UnlimitedFastNBLVerletList
{
  public:
    using traits_type         = traitsT;
    using system_type         = System<traits_type>;
    using boundary_type       = typename traits_type::boundary_type;
    using real_type           = typename traits_type::real_type;
    using coordinate_type     = typename traits_type::coordinate_type;
    using exclusion_list_type = ExclusionList;

    using parameter_type      = parameterT;
    using neighbor_list_type  = NeighborList<parameter_type>;
    using neighbor_type       = typename neighbor_list_type::neighbor_type;
    using range_type          = typename neighbor_list_type::range_type;

  public:
    UnlimitedFastNBLVerletList() : margin_(0.5), current_margin_(-1.0){}
    explicit UnlimitedFastNBLVerletList(const real_type mgn)
        : margin_(mgn), current_margin_(-1.0)
    {}

    ~UnlimitedFastNBLVerletList() = default;
    UnlimitedFastNBLVerletList(UnlimitedFastNBLVerletList const&) = default;
    UnlimitedFastNBLVerletList(UnlimitedFastNBLVerletList &&)     = default;
    UnlimitedFastNBLVerletList& operator=(UnlimitedFastNBLVerletList const&) = default;
    UnlimitedFastNBLVerletList& operator=(UnlimitedFastNBLVerletList &&)     = default;

    bool valid() const noexcept
    {
        return current_margin_ >= 0.0;
    }

    template<typename PotentialT>
    void initialize(const system_type& sys, const PotentialT& pot)
    {
        this->set_cutoff(pot.max_cutoff_length());
        this->exclusion_.make(sys, pot);
        this->make(sys, pot);
        return;
    }

    template<typename PotentialT>
    void make(const system_type& sys, const PotentialT& pot)
    {
        static_assert(std::is_same<
            typename PotentialT::parameter_type, parameter_type>::value, "");

        this->neighbors_.clear();

        // `participants` is a list that contains indices of particles that are
        // related to the potential.
        const auto& participants = pot.participants();

        const real_type rc  = cutoff_ * (1. + margin_);
        const real_type rc2 = rc * rc;
        // set r_c = 3*L
        //              3            2           1
        //            2109 8765 4321 0987 6543 2109 8765 4321
        // mask_0 = 0b0010'0000'0000'1000'0000'0010'0000'0000
        //        =      2    0    0    8    0    2    0    0
        // mask_1 = 0b0001'1111'1100'0111'1111'0001'1111'1100
        //        =      1    F    C    7    F    1    F    C
        constexpr std::uint32_t mask_0  = 0x20080200;
        constexpr std::uint32_t mask_1  = 0x1FC7F1FC;
        constexpr std::uint32_t mask_9  = 0x000001FF; // extract last 9 bits
        constexpr real_type     ratio = 1.0 / 3.0;

        const real_type L  = rc * ratio;
        const real_type rL = real_type(1.0) / L; // reciprocal L

        this->lattice_coords_.resize(participants.size());
        for(std::size_t i=0; i<participants.size(); ++i)
        {
            const auto r = sys.position(participants[i]) * rL;
            this->lattice_coords_[i] =
                (std::uint32_t(std::floor(math::Z(r))) & mask_9 << 20) +
                (std::uint32_t(std::floor(math::Y(r))) & mask_9 << 10) +
                (std::uint32_t(std::floor(math::X(r))) & mask_9);
        }

        std::vector<neighbor_type> partner;
        for(std::size_t idx=0; idx<participants.size(); ++idx)
        {
            partner.clear();
            const auto  i  = participants[idx];
            const auto  ri = this->lattice_coords_[idx];
            const auto& qi = sys.position(i);

            for(std::size_t jdx=idx+1; jdx<participants.size(); ++jdx)
            {
                const auto j = participants[jdx];
                if(this->exclusion_.is_excluded(i, j))
                {
                    continue;
                }
                const auto rj = this->lattice_coords_[idx];

                const std::uint32_t D0 = (ri + mask_0 - rj) & mask_1;
                const std::uint32_t D1 = (rj + mask_0 - ri) & mask_1;

                const std::uint32_t R0 =  D0       & D1;
                const std::uint32_t R1 = (D0 << 1) & D1;
                const std::uint32_t R2 = (D0 >> 1) & D1;

                if(R0 == 0 && R1 == 0 && R2 == 0)
                {
                    const auto& qj = sys.position(j);
                    if(math::length_sq(sys.adjust_direction(qj - qi)) < rc2)
                    {
                        partner.emplace_back(j, pot.prepare_params(i, j));
                    }
                }
            }
            // because j is searched sequencially, sorting is not needed.
            this->neighbors_.add_list_for(i, partner.begin(), partner.end());
        }
        this->current_margin_ = cutoff_ * margin_;
        return ;
    }

    template<typename PotentialT>
    void update(const real_type dmargin, const system_type& sys,
                const PotentialT& pot)
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(sys, pot);
        }
        return ;
    }

    real_type cutoff() const noexcept {return this->cutoff_;}
    real_type margin() const noexcept {return this->margin_;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    void set_cutoff(const real_type c) noexcept {this->cutoff_ = c;}
    void set_margin(const real_type m) noexcept {this->margin_ = m;}

  private:

    real_type      cutoff_;
    real_type      margin_;
    real_type      current_margin_;

    std::vector<std::uint32_t> lattice_coords_;
    exclusion_list_type        exclusion_;
    neighbor_list_type         neighbors_;
};
} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
