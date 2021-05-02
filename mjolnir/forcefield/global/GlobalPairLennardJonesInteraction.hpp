#ifndef MJOLNIR_INTEARACTION_GLOBAL_PAIR_LENNARD_JONES_INTEARACTION_HPP
#define MJOLNIR_INTEARACTION_GLOBAL_PAIR_LENNARD_JONES_INTEARACTION_HPP
#include <mjolnir/forcefield/global/GlobalPairInteraction.hpp>
#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <memory>

namespace mjolnir
{

// It is a specialization of GlobalPairInteraction for LennardJonesPotential.
// In the case of LennardJonesPotential, we can omit `sqrt` call that is
// normally used to calculate distance because we only needs the squared distance.
template<typename realT, template<typename, typename> class boundaryT>
class GlobalPairInteraction<
    SimulatorTraits<realT, boundaryT>,
    LennardJonesPotential<realT>
    > final : public GlobalInteractionBase<SimulatorTraits<realT, boundaryT>>
{
  public:

    using traits_type         = SimulatorTraits<realT, boundaryT>;
    using base_type           = GlobalInteractionBase<traits_type>;
    using real_type           = typename base_type::real_type;
    using coordinate_type     = typename base_type::coordinate_type;
    using system_type         = typename base_type::system_type;
    using topology_type       = typename base_type::topology_type;
    using boundary_type       = typename base_type::boundary_type;
    using potential_type      = LennardJonesPotential<real_type>;
    using partition_type      = SpatialPartition<traits_type, potential_type>;
    using parameter_list_type = GlobalParameterList<traits_type, potential_type>;

  public:

    GlobalPairInteraction(parameter_list_type&& para, partition_type&& part)
        : parameters_(std::move(para)), partition_(std::move(part))
    {}
    ~GlobalPairInteraction() override {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->parameters_.initialize(sys, topol);
        this->partition_ .initialize(sys, parameters_);
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->parameters_.update(sys, topol);
        this->partition_.initialize(sys, this->parameters_);
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.reduce_margin(dmargin, sys, this->parameters_);
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->partition_.scale_margin(scale, sys, this->parameters_);
        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        const auto cutoff_ratio    = potential_type::cutoff_ratio;
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;

        const auto leading_participants = this->parameters_.leading_participants();
        for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
        {
            const auto i = leading_participants[idx];

            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j   = ptnr.index;
                const auto& pot = ptnr.potential();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = pot.sigma() * pot.sigma();
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type epsilon = pot.epsilon();

                const real_type rcp_l_sq = 1 / l_sq;
                const real_type s2l2 = sigma_sq * rcp_l_sq;
                const real_type s6l6 = s2l2 * s2l2 * s2l2;

                const coordinate_type f = rij *
                    (24 * epsilon * (s6l6 - 2 * s6l6 * s6l6) * rcp_l_sq);

                sys.force(i) += f;
                sys.force(j) -= f;
            }
        }
        return ;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E(0);

        const auto cutoff_ratio    = potential_type::cutoff_ratio;
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;
        const auto coef_at_cutoff  = potential_type::coef_at_cutoff;

        const auto leading_participants = this->parameters_.leading_participants();
        for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
        {
            const auto i = leading_participants[idx];

            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j   = ptnr.index;
                const auto& pot = ptnr.potential();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = pot.sigma() * pot.sigma();
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type epsilon = pot.epsilon();

                const real_type s2l2 = sigma_sq / l_sq;
                const real_type s6l6 = s2l2 * s2l2 * s2l2;

                E += 4 * epsilon * (s6l6 * s6l6 - s6l6 - coef_at_cutoff);
            }
        }
        return E;
    }

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        const auto cutoff_ratio    = potential_type::cutoff_ratio;
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;
        const auto coef_at_cutoff  = potential_type::coef_at_cutoff;

        real_type energy = 0;
        const auto leading_participants = this->parameters_.leading_participants();
        for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
        {
            const auto i = leading_participants[idx];

            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j   = ptnr.index;
                const auto& pot = ptnr.potential();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = pot.sigma() * pot.sigma();
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type epsilon = pot.epsilon();

                const real_type rcp_l_sq = 1 / l_sq;
                const real_type s2l2 = sigma_sq * rcp_l_sq;
                const real_type s6l6 = s2l2 * s2l2 * s2l2;

                energy += 4 * epsilon * (s6l6 * s6l6 - s6l6 - coef_at_cutoff);

                const coordinate_type f = rij *
                    (24 * epsilon * (s6l6 - 2 * s6l6 * s6l6) * rcp_l_sq);

                sys.force(i) += f;
                sys.force(j) -= f;
            }
        }
        return energy;
    }

    std::string name() const override {return "GlobalPairLennardJones";}

    base_type* clone() const override
    {
        return new GlobalPairInteraction(
                parameter_list_type(parameters_), partition_type(partition_));
    }

  private:

    parameter_list_type parameters_;
    partition_type partition_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize BondAngleInteraction with LocalPotentials
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
// L-J
extern template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float >>;
extern template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
