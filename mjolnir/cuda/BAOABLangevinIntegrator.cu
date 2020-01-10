#include <mjolnir/cuda/BAOABLangevinIntegrator.hpp>
#include <mjolnir/cuda/Vector.hpp>
#include <thrust/tuple.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <cuda.h>
#include <curand.h>

namespace mjolnir
{

// ---------------------------------------------------------------------------
// To reduce implementation cost, first implement template function.
// This cannot be written in a header file because it uses cuda device
// operations and the header file would be included from .cpp file that will be
// compiled by gcc. When this template function is instanciated via cpp file,
// gcc cannot find the way to compile cuda operations and stop with an error.
// Therefore, we need to write all the definitions in a .cu file.

namespace detail
{

template<typename realT, template<typename, typename> class boundaryT>
realT cuda_BAOABLangevinIntegrator_step1(
    BAOABLangevinIntegrator<CUDASimulatorTraits<realT, boundaryT>>& integ,
    System                 <CUDASimulatorTraits<realT, boundaryT>>& sys,
    RandomNumberGenerator  <CUDASimulatorTraits<realT, boundaryT>>& rng) noexcept
{
    using namespace ::mjolnir::math; // for real4 operators
    using traits_type     = CUDASimulatorTraits<realT, UnlimitedBoundary>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    rng.fill_with_gaussian(integ.random_forces_device());
    thrust::device_ptr<const real_type> rnd = integ.random_forces_device().data();

    const real_type kB     = physics::constants<real_type>::kB();
    const real_type T_ref  = sys.attribute("temperature");
    const real_type kBT    = kB * T_ref;
    const real_type halfdt = integ.halfdt();

    const thrust::device_ptr<const real_type> expgts  = integ.exp_gamma_dt_device().data();
    const thrust::device_ptr<const real_type> noises  = integ.noise_coeff_device().data();
    const thrust::device_ptr<const real_type> rmasses = sys.rmasses_device().data();
    const thrust::device_ptr<coordinate_type> ps      = sys.positions_device().data();
    const thrust::device_ptr<coordinate_type> vs      = sys.velocities_device().data();
    const thrust::device_ptr<coordinate_type> fs      = sys.forces_device().data();

    const auto boundary = sys.boundary();

    return thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<std::size_t>(0),
        thrust::counting_iterator<std::size_t>(sys.size()),
        [halfdt, kBT, rnd, expgts, noises, rmasses, ps, vs, fs, boundary]
        __device__ (const std::size_t idx)
        {
            const std::size_t offset = idx * 3;
            coordinate_type R;
            R.x = rnd.get()[offset+0];
            R.y = rnd.get()[offset+1];
            R.z = rnd.get()[offset+2];
            R.w = 0;

            const real_type  expgt = expgts .get()[idx];
            const real_type  noise = noises .get()[idx];
            const real_type  rm    = rmasses.get()[idx];
            coordinate_type& p     = ps.get()[idx];
            coordinate_type& v     = vs.get()[idx];
            coordinate_type& f     = fs.get()[idx];
            coordinate_type  dp{0,0,0,0};

            v  += halfdt * rm * f; // calc v(n+1/3)
            dp += halfdt * v;      // calc p(n+1/2)
            v  *= expgt;
            v  += noise  * R;      // calc v(n+2/3)
            dp += halfdt * v;      // calc p(n+1)

            // update p(n) -> p(n+1);
            p  = boundary.adjust_position(p + dp);

            // reset force
            f = coordinate_type{0,0,0,0};

            // collect largest displacement
            return dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;

        }, 0.0, thrust::maximum<double>());
}

template<typename traitsT>
void cuda_BAOABLangevinIntegrator_step2(BAOABLangevinIntegrator<traitsT>& integ,
                                        System<traitsT>& sys) noexcept
{
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    const real_type halfdt = integ.halfdt();
    const thrust::device_ptr<const real_type> rmasses = sys.rmasses_device().data();
    const thrust::device_ptr<coordinate_type> vs      = sys.velocities_device().data();
    const thrust::device_ptr<coordinate_type> fs      = sys.forces_device().data();

    // calc v(n+2/3) -> v(n+1)
    thrust::for_each(thrust::device,
        thrust::counting_iterator<std::size_t>(0),
        thrust::counting_iterator<std::size_t>(sys.size()),
        [halfdt, rmasses, vs, fs] __device__ (const std::size_t idx)
        {
            vs.get()[idx] += halfdt * rmasses[idx] * fs[idx];
            return;
        });
    return ;
}

template<typename traitsT>
void cuda_BAOABLangevinIntegrator_update_parameters(
        BAOABLangevinIntegrator<traitsT>& integ, const System<traitsT>& sys) noexcept
{
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    const real_type dt  = integ.delta_t();
    const real_type kBT = physics::constants<real_type>::kB() * integ.temperature();

    const thrust::device_ptr<const real_type> rmasses = sys.rmasses_device().data();
    const thrust::device_ptr<const real_type> gammas  = integ.gammas_device().data();
    const thrust::device_ptr<real_type>       expgts  = integ.exp_gamma_dt_device().data();
    const thrust::device_ptr<real_type>       noises  = integ.noise_coeff_device().data();

    thrust::for_each(thrust::device,
        thrust::counting_iterator<std::size_t>(0),
        thrust::counting_iterator<std::size_t>(sys.size()),
        [dt, kBT, rmasses, gammas, expgts, noises]
        __device__ (const std::size_t idx)
        {
            const real_type rmass    = rmasses[idx];
            const real_type gamma_dt = -real_type(1) * gammas[idx] * dt;
            expgts[idx] = ::exp(gamma_dt);
            noises[idx] = ::sqrt(kBT * (real_type(1) - ::exp(real_type(2) * gamma_dt)) * rmass);
            return ;
        });
    return ;
}
} // detail

// ----------------------------------------------------------------------------
// definitions of BAOABLangevinIntegrator::step()

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    const real_type largest_disp2 =
        detail::cuda_BAOABLangevinIntegrator_step1(*this, sys, rng);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff.reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(p(n+1))
    ff.calc_force(sys);

    // calc v(n+2/3) -> v(n+1)
    detail::cuda_BAOABLangevinIntegrator_step2(*this, sys);

    // XXX Note that the updated configurations are not pulled here.
    //     Almost all the computations are done on GPU in this implemnetation.
    //     Observers require the values, so they pull it back to CPU.

    return time + dt_;
}

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    const real_type largest_disp2 =
        detail::cuda_BAOABLangevinIntegrator_step1(*this, sys, rng);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff.reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(p(n+1))
    ff.calc_force(sys);

    // calc v(n+2/3) -> v(n+1)
    detail::cuda_BAOABLangevinIntegrator_step2(*this, sys);

    // XXX Note that the updated configurations are not pulled here.
    //     Almost all the computations are done on GPU in this implemnetation.
    //     Observers require the values, so they pull it back to CPU.

    return time + dt_;
}

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    const real_type largest_disp2 =
        detail::cuda_BAOABLangevinIntegrator_step1(*this, sys, rng);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff.reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(p(n+1))
    ff.calc_force(sys);

    // calc v(n+2/3) -> v(n+1)
    detail::cuda_BAOABLangevinIntegrator_step2(*this, sys);

    // XXX Note that the updated configurations are not pulled here.
    //     Almost all the computations are done on GPU in this implemnetation.
    //     Observers require the values, so they pull it back to CPU.

    return time + dt_;
}

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    const real_type largest_disp2 =
        detail::cuda_BAOABLangevinIntegrator_step1(*this, sys, rng);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff.reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(p(n+1))
    ff.calc_force(sys);

    // calc v(n+2/3) -> v(n+1)
    detail::cuda_BAOABLangevinIntegrator_step2(*this, sys);

    // XXX Note that the updated configurations are not pulled here.
    //     Almost all the computations are done on GPU in this implemnetation.
    //     Observers require the values, so they pull it back to CPU.

    return time + dt_;
}

// ---------------------------------------------------------------------------
// definitions of BAOABLangevinIntegrator<CUDA>::initialize();

template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&)
{
    // calculate parameters for each particles
    this->update(sys);

    // zero-clear force
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }
    // also on the device.
    sys.forces_device() = sys.forces_host();

    // calculate the current force.
    ff.calc_force(sys);
    return;
}
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&)
{
    // calculate parameters for each particles
    this->update(sys);

    // zero-clear force
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }
    // also on the device.
    sys.forces_device() = sys.forces_host();

    // calculate the current force.
    ff.calc_force(sys);
    return;
}
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&)
{
    // calculate parameters for each particles
    this->update(sys);

    // zero-clear force
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }
    // also on the device.
    sys.forces_device() = sys.forces_host();

    // calculate the current force.
    ff.calc_force(sys);
    return;
}
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&)
{
    // calculate parameters for each particles
    this->update(sys);

    // zero-clear force
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }
    // also on the device.
    sys.forces_device() = sys.forces_host();

    // calculate the current force.
    ff.calc_force(sys);
    return;
}

// ---------------------------------------------------------------------------
// definitions of BAOABLangevinIntegrator<CUDA>::reset_parameters();

template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>
    >::reset_parameters(const system_type& sys)
{
    detail::cuda_BAOABLangevinIntegrator_update_parameters(*this, sys);
    return;
}

template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>
    >::reset_parameters(const system_type& sys)
{
    detail::cuda_BAOABLangevinIntegrator_update_parameters(*this, sys);
    return;
}

template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>
    >::reset_parameters(const system_type& sys)
{
    detail::cuda_BAOABLangevinIntegrator_update_parameters(*this, sys);
    return;
}

template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>
    >::reset_parameters(const system_type& sys)
{
    detail::cuda_BAOABLangevinIntegrator_update_parameters(*this, sys);
    return;
}

// ----------------------------------------------------------------------------
// any constructors and `operator=`s are the same as the one that will be
// generated by `= default`, but it should be here because those uses
// `thrust::device_vector`s copy/move constructors and `operator=`s that uses
// device instructions.
//     Since all the implementations are the same, it uses macro.

#define GENERATE_BAOAB_LANGEVIN_INTEGRATOR_CONSTRUCTOR_INSTANCES(REAL, BOUNDARY)          \
    template<>                                                                            \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>::BAOABLangevinIntegrator(\
            const real_type dt, const std::vector<real_type>& gamma)                      \
        : dt_(dt), halfdt_(dt / 2), gammas_h_(gamma.begin(), gamma.end()),                \
          gammas_(gammas_h_), exp_gamma_dt_(gammas_.size()),                              \
          noise_coeff_ (gammas_.size()), random_forces_(gammas_.size() * 3)               \
    {} /**/

GENERATE_BAOAB_LANGEVIN_INTEGRATOR_CONSTRUCTOR_INSTANCES(double, UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_CONSTRUCTOR_INSTANCES(float,  UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_CONSTRUCTOR_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_CONSTRUCTOR_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_BAOAB_LANGEVIN_INTEGRATOR_CONSTRUCTOR_INSTANCES

#define GENERATE_BAOAB_LANGEVIN_INTEGRATOR_DESTRUCTOR_INSTANCES(REAL, BOUNDARY)             \
    template<>                                                                              \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>::~BAOABLangevinIntegrator()\
    {} /**/

GENERATE_BAOAB_LANGEVIN_INTEGRATOR_DESTRUCTOR_INSTANCES(double, UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_DESTRUCTOR_INSTANCES(float,  UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_DESTRUCTOR_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_DESTRUCTOR_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_BAOAB_LANGEVIN_INTEGRATOR_DESTRUCTOR_INSTANCES

#define GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_CONSTRUCTOR_INSTANCES(REAL, BOUNDARY)     \
    template<>                                                                            \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>::BAOABLangevinIntegrator(\
            const BAOABLangevinIntegrator& other)                                         \
        : dt_(other.dt_), halfdt_(other.halfdt_), gammas_h_(other.gammas_h_),             \
          gammas_(other.gammas_), exp_gamma_dt_(other.exp_gamma_dt_),                     \
          noise_coeff_ (other.noise_coeff_), random_forces_(other.random_forces_)         \
    {} /**/

GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_CONSTRUCTOR_INSTANCES(double, UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_CONSTRUCTOR_INSTANCES(float,  UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_CONSTRUCTOR_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_CONSTRUCTOR_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_CONSTRUCTOR_INSTANCES

#define GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_CONSTRUCTOR_INSTANCES(REAL, BOUNDARY)      \
    template<>                                                                             \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>::BAOABLangevinIntegrator( \
            BAOABLangevinIntegrator&& other)                                               \
        : dt_(other.dt_), halfdt_(other.halfdt_), gammas_h_(std::move(other.gammas_h_)),   \
          gammas_(std::move(other.gammas_)), exp_gamma_dt_(std::move(other.exp_gamma_dt_)),\
          noise_coeff_ (std::move(other.noise_coeff_)),                                    \
          random_forces_(std::move(other.random_forces_))                                  \
    {} /**/

GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_CONSTRUCTOR_INSTANCES(double, UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_CONSTRUCTOR_INSTANCES(float,  UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_CONSTRUCTOR_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_CONSTRUCTOR_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_CONSTRUCTOR_INSTANCES

#define GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_ASSIGNMENT_INSTANCES(REAL, BOUNDARY)\
    template<>                                                                      \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>&                   \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>::operator=(        \
            const BAOABLangevinIntegrator& other)                                   \
    {                                                                               \
        this->dt_            = other.dt_;                                           \
        this->halfdt_        = other.halfdt_;                                       \
        this->gammas_h_      = other.gammas_h_;                                     \
        this->gammas_        = other.gammas_;                                       \
        this->exp_gamma_dt_  = other.exp_gamma_dt_;                                 \
        this->noise_coeff_   = other.noise_coeff_;                                  \
        this->random_forces_ = other.random_forces_;                                \
        return *this;                                                               \
    } /**/

GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_ASSIGNMENT_INSTANCES(double, UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_ASSIGNMENT_INSTANCES(float,  UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_ASSIGNMENT_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_ASSIGNMENT_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_BAOAB_LANGEVIN_INTEGRATOR_COPY_ASSIGNMENT_INSTANCES

#define GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_ASSIGNMENT_INSTANCES(REAL, BOUNDARY)\
    template<>                                                                      \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>&                   \
    BAOABLangevinIntegrator<CUDASimulatorTraits<REAL, BOUNDARY>>::operator=(        \
            BAOABLangevinIntegrator&& other)                                        \
    {                                                                               \
        this->dt_            = std::move(other.dt_);                                \
        this->halfdt_        = std::move(other.halfdt_);                            \
        this->gammas_h_      = std::move(other.gammas_h_);                          \
        this->gammas_        = std::move(other.gammas_);                            \
        this->exp_gamma_dt_  = std::move(other.exp_gamma_dt_);                      \
        this->noise_coeff_   = std::move(other.noise_coeff_);                       \
        this->random_forces_ = std::move(other.random_forces_);                     \
        return *this;                                                               \
    } /**/

GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_ASSIGNMENT_INSTANCES(double, UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_ASSIGNMENT_INSTANCES(float,  UnlimitedBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_ASSIGNMENT_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_ASSIGNMENT_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_BAOAB_LANGEVIN_INTEGRATOR_MOVE_ASSIGNMENT_INSTANCES

// ---------------------------------------------------------------------------
// class template instanciation

template class BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
