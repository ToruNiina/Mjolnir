#ifndef MJOLNIR_CORE_GJ_F_NPT_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CORE_GJ_F_NPT_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// "A simple and effective Verlet-type algorithm for simulating Langevin dynamics"
// Niels Gronbech-Jensen and Oded Farago, Molecular Physics, 2013, Vol.111 No.8, 983-991

template<typename traitsT>
class GJFNPTLangevinIntegrator
{
  public:
    using traits_type     = traitsT;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traitsT>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traitsT>>;
    using rng_type        = RandomNumberGenerator<traits_type>;
    using remover_type    = SystemMotionRemover<traits_type>;

  public:

    GJFNPTLangevinIntegrator(const real_type dt,
            const real_type box_alpha, const real_type box_Q,
            const real_type volume_velo,
            std::vector<real_type>&& alpha, remover_type&& remover)
        : dt_(dt),
          box_alpha_(box_alpha), box_Q_(box_Q), volume_velocity_(volume_velo),
          alphas_(std::move(alpha)), betas_(alphas_.size()), bs_(alphas_.size()),
          remover_(std::move(remover))
    {}
    ~GJFNPTLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type&)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if( ! is_cuboidal_periodic_boundary<boundary_type>::value ||
            math::X(sys.boundary().lower_bound()) != real_type(0) ||
            math::Y(sys.boundary().lower_bound()) != real_type(0) ||
            math::Z(sys.boundary().lower_bound()) != real_type(0))
        {
            MJOLNIR_LOG_ERROR("NPT Integrator requires periodic boundary of which"
                              " lower bound is at the origin (0, 0, 0).");
            throw_exception<std::runtime_error>("NPT Integrator requires periodic"
                    " boundary of which lower bound is at the origin (0, 0, 0).");
        }

        if( ! ff->constraint().empty())
        {
            MJOLNIR_LOG_WARN("G-JF NPT langevin integrator does not support"
                "  constraints. [[forcefields.constraint]] will be ignored.");
        }

        // calculate parameters for each particles
        this->update(sys);

        // if loaded from MsgPack, we can skip it.
        if( ! sys.force_initialized())
        {
            // calculate force
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
            }
            ff->calc_force(sys);
        }

        this->virial_ = 0;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->virial_ += math::dot_product(sys.force(i), sys.position(i));
        }

        assert(math::X(sys.boundary().lower_bound()) == 0);
        assert(math::Y(sys.boundary().lower_bound()) == 0);
        assert(math::Z(sys.boundary().lower_bound()) == 0);

        assert(math::X(sys.boundary().upper_bound()) > 0);
        assert(math::Y(sys.boundary().upper_bound()) > 0);
        assert(math::Z(sys.boundary().upper_bound()) > 0);

        const auto width = sys.boundary().upper_bound();
        volume_ = math::X(width) * math::Y(width) * math::Z(width);

        const auto kBT = physics::constants<real_type>::kB() * this->temperature_;

        const auto internal_pressure = (virial_ / 3 + sys.size() * kBT) / volume_;
        sys.attribute("internal_pressure") = internal_pressure;

        this->volume_force_ = internal_pressure - pressure_;

        return;
    }

    real_type step(const real_type, system_type&, forcefield_type&, rng_type&);

    void update(const system_type& sys)
    {
        if(!sys.has_attribute("temperature"))
        {
            throw std::out_of_range("mjolnir::GJFNPTLangevinIntegrator: "
                "Langevin Integrator requires reference temperature, but "
                "`temperature` is not found in `system.attribute`.");
        }
        this->temperature_ = sys.attribute("temperature");
        this->reset_parameters(sys);
        return;
    }

    real_type delta_t() const noexcept {return dt_;}
    std::vector<real_type> const& parameters() const noexcept {return alphas_;}

  private:

    void reset_parameters(const system_type& sys) noexcept
    {
        const auto kBT = physics::constants<real_type>::kB() * this->temperature_;

        box_1_over_2Q_ = real_type(0.5) / box_Q_;

        box_beta_ = std::sqrt(real_type(2) * box_alpha_ * kBT * dt_);
        box_a_    = (real_type(1) - box_alpha_ * dt_ * box_1_over_2Q_) /
                    (real_type(1) + box_alpha_ * dt_ * box_1_over_2Q_);
        box_b_    = real_type(1) /
                    (real_type(1) + box_alpha_ * dt_ * box_1_over_2Q_);

        betas_    .resize(sys.size());
        bs_       .resize(sys.size());
        as_       .resize(sys.size());
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto alpha   = this->alphas_.at(i);
            const auto m       = sys.mass(i);
            this->as_.at(i)    = (1.0 - alpha * dt_ * 0.5 / m) /
                                 (1.0 + alpha * dt_ * 0.5 / m);
            this->bs_.at(i)    = 1.0 / (1.0 + alpha * dt_ * 0.5 / m);
            this->betas_.at(i) = std::sqrt(2 * alpha * kBT * dt_);
        }
        return;
    }

    coordinate_type gen_R(rng_type& rng) noexcept
    {
        const auto x = rng.gaussian();
        const auto y = rng.gaussian();
        const auto z = rng.gaussian();
        return math::make_coordinate<coordinate_type>(x, y, z);
    }

  private:
    real_type dt_;
    real_type temperature_; // reference T
    real_type pressure_;    // reference V
    real_type virial_;      // just sum f*r

    // for the notation, see the paper.

    real_type box_alpha_;
    real_type box_beta_;
    real_type box_a_;
    real_type box_b_;
    real_type box_Q_;
    real_type box_1_over_2Q_;

    real_type volume_;
    real_type volume_velocity_;
    real_type volume_force_;

    std::vector<real_type> alphas_;
    std::vector<real_type> betas_;
    std::vector<real_type> as_;
    std::vector<real_type> bs_;

    remover_type remover_;
};

template<typename traitsT>
typename GJFNPTLangevinIntegrator<traitsT>::real_type
GJFNPTLangevinIntegrator<traitsT>::step(const real_type time,
        system_type& sys, forcefield_type& ff, rng_type& rng)
{
    // update volume^(n+1)
    const auto prev_volume = this->volume_; // volume^(n)

    this->volume_ += box_b_ * dt_ * volume_velocity_ +
                     box_b_ * dt_ * dt_ * box_1_over_2Q_ * volume_force_ +
                     box_b_ * dt_ * box_1_over_2Q_ * box_beta_ * rng.gaussian();

    // update position (and velocity, partially)

    const real_type l_prev = std::cbrt(prev_volume);
    const real_type l_curr = std::cbrt(volume_);
    const real_type l_scale1 = l_curr / l_prev;
    const real_type l_scale2 = 2 * l_curr / (l_curr + l_prev);

    sys.boundary().set_upper_bound(sys.boundary().upper_bound() * l_scale1);

    real_type largest_disp2(0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto rm = sys.rmass(i);  // reciprocal mass
        auto&       p = sys.position(i);
        auto&       v = sys.velocity(i);
        auto&       f = sys.force(i);
        const auto& b = this->bs_[i];
        const auto& a = this->as_[i];

        const auto beta = this->gen_R(rng) * betas_[i]; // gen beta^(n+1)

        const auto dp = l_scale2 * b * dt_ *
            (v + (dt_ * rm * real_type(0.5)) * f + (rm * real_type(0.5)) * beta);

        p *= l_scale1; // this moves p into the boundary (box has been changed)
        p  = sys.adjust_position(p + dp); // move particle

        // update velocity partially
        v *= a;
        v += (a * dt_ * rm * real_type(0.5)) * f + (b * rm) * beta;

        // reset force
        f = math::make_coordinate<coordinate_type>(0, 0, 0);

        // collect largest displacement
        largest_disp2 = std::max(largest_disp2, math::length_sq(dp));
    }

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff->scale_margin(l_scale1, sys);
    ff->reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc force^(n+1)

    ff->calc_force(sys);

    // calc velocity^(n+1)

    this->virial_ = 0;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto  rm   = sys.rmass(i);  // reciprocal math
        const auto& p    = sys.position(i);
        auto&       v    = sys.velocity(i);
        const auto& f    = sys.force(i);

        v += (dt_ * rm * real_type(0.5)) * f;

        this->virial_ += math::dot_product(f, p);
    }

    // update volume force/velocity ^(n+1)

    const auto kBT = physics::constants<real_type>::kB() * this->temperature_;

    this->volume_force_ = (virial_ / 3 + sys.size() * kBT) / volume_ - pressure_;
    this->volume_velocity_ += dt_ * box_1_over_2Q_ * volume_force_;

    sys.attribute("internal_pressure") = (virial_ / 3 + sys.size() * kBT) / volume_;

    // remove net rotation/translation
    remover_.remove(sys);

    return time + dt_;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GJFNPTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GJFNPTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif// MJOLNIR_CORE_GJ_F_NPT_LANGEVIN_INTEGRATOR_HPP
