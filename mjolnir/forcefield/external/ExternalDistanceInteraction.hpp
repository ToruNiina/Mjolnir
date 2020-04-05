#ifndef MJOLNIR_INTERACTION_EXTERNAL_DISTANCE_INTERACTION_HPP
#define MJOLNIR_INTERACTION_EXTERNAL_DISTANCE_INTERACTION_HPP
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>

namespace mjolnir
{

/*! @brief Interaction between particle and a region based on their distance. *
 *  @details shapeT represents the shape. It provides a method to calculate   *
 *           distance between particle and the shape, force direction, and    *
 *           neighbor-list.                                                   */
template<typename traitsT, typename potentialT, typename shapeT>
class ExternalDistanceInteraction final
    : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using potential_type  = potentialT;
    using shape_type      = shapeT;
    using base_type       = ExternalForceInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

  public:

    ExternalDistanceInteraction(shape_type&& shape, potential_type&& pot)
        : shape_(std::move(shape)), potential_(std::move(pot))
    {}
    ~ExternalDistanceInteraction() override {}

    // calculate force, update spatial partition (reduce margin) inside.
    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(system_type const&) const noexcept override;

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        this->potential_.update(sys); // update system parameters
        this->shape_.initialize(sys, this->potential_);
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys) override
    {
        this->potential_.update(sys); // update system parameters
        this->shape_.initialize(sys, this->potential_);
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->shape_.reduce_margin(dmargin, sys);
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->shape_.scale_margin(scale, sys);
    }

    std::string name() const override
    {return "ExternalDistance:"_s + potential_.name();}

    base_type* clone() const override
    {
        return new ExternalDistanceInteraction(
                shape_type(shape_), potential_type(potential_));
    }

  private:

    shape_type     shape_;
    potential_type potential_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potT, typename spaceT>
void ExternalDistanceInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys) const noexcept
{
    for(std::size_t i : this->shape_.neighbors())
    {
        const auto& ri = sys.position(i);

        const real_type dist = this->shape_.calc_distance(ri, sys.boundary());
        const real_type dV   = this->potential_.derivative(i, dist);
        if(dV == 0.0){continue;}

        const auto f = shape_.calc_force_direction(ri, sys.boundary());
        sys.force(i) += -dV * f;
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename ExternalDistanceInteraction<traitsT, potT, spaceT>::real_type
ExternalDistanceInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(std::size_t i : this->shape_.neighbors())
    {
        const auto&    ri = sys.position(i);
        const real_type d = this->shape_.calc_distance(ri, sys.boundary());
        E += this->potential_.potential(i, d);
    }
    return E;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/forcefield/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/forcefield/external/ImplicitMembranePotential.hpp>
#include <mjolnir/forcefield/external/LennardJonesWallPotential.hpp>

namespace mjolnir
{

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ExcludedVolumeWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ExcludedVolumeWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;


extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, LennardJonesWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, LennardJonesWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;


extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, PositiveXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , PositiveXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, NegativeXAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , NegativeXAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, PositiveYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , PositiveYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, NegativeYAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , NegativeYAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, PositiveZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , PositiveZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ExternalDistanceInteraction<SimulatorTraits<double, UnlimitedBoundary       >, ImplicitMembranePotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, ImplicitMembranePotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  UnlimitedBoundary       >>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ImplicitMembranePotential<double>, NegativeZAxisAlignedPlane<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ExternalDistanceInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ImplicitMembranePotential<float> , NegativeZAxisAlignedPlane<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;


} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif//MJOLNIR_BOX_INTEARACTION_BASE
