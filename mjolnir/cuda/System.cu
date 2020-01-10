#include <mjolnir/cuda/System.hpp>
#include <mjolnir/cuda/Vector.hpp>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_ptr.h>
#include <cuda.h>

namespace mjolnir
{

// ---------------------------------------------------------------------------
// definitions of BAOABLangevinIntegrator<CUDA>::step();

namespace detail
{
template<typename traitsT>
void cuda_system_initialize(
        System<traitsT>& sys, RandomNumberGenerator<traitsT>& rng)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const real_type kB    = physics::constants<real_type>::kB();
    const real_type T_ref = sys.attribute("temperature");
    const real_type kBT   = kB * T_ref;

    MJOLNIR_LOG_NOTICE("generating velocity with T = ", T_ref, "...");

    // generate gaussians using cuRAND

    thrust::device_vector<real_type> gaussians(sys.size() * 3);
    rng.fill_with_gaussian(gaussians);

    thrust::device_ptr<const real_type> rmass = sys.rmasses_device().data();
    thrust::device_ptr<coordinate_type> velo  = sys.velocities_device().data();
    thrust::device_ptr<const real_type> rnd   = gaussians.data();

    // generate Maxwell-Boltzmann distribution
    thrust::for_each(thrust::device,
        thrust::counting_iterator<std::size_t>(0),
        thrust::counting_iterator<std::size_t>(sys.size()),
        [kBT, rmass, velo, rnd] __device__ (const std::size_t idx)
        {
            const std::size_t offset = idx * 3;
            const auto R = math::make_coordinate<coordinate_type>(
                    rnd[offset], rnd[offset+1], rnd[offset+2]);

            velo[idx] = ::sqrt(kBT * rmass[idx]) * R;
            return;
        });

    // pull back generated velocities to system
    sys.velocities_host() = sys.velocities_device();

    MJOLNIR_LOG_NOTICE("done.");
    return;
}

} // detail

template<>
void System<CUDASimulatorTraits<double, UnlimitedBoundary>>::initialize(rng_type& rng)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // copy initial conditions from CPU to GPU
    this->masses_d_     = this->masses_    ;
    this->rmasses_d_    = this->rmasses_   ;
    this->positions_d_  = this->positions_ ;
    this->velocities_d_ = this->velocities_;
    this->forces_d_     = this->forces_    ;

    if(this->velocity_initialized_)
    {
        MJOLNIR_LOG_NOTICE(
            "velocity is already given, nothing to initialize in System");
        return ;
    }
    assert(this->has_attribute("temperature"));

    detail::cuda_system_initialize(*this, rng);
    return;
}

template<>
void System<CUDASimulatorTraits<float, UnlimitedBoundary>>::initialize(rng_type& rng)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // copy initial conditions from CPU to GPU
    this->masses_d_     = this->masses_    ;
    this->rmasses_d_    = this->rmasses_   ;
    this->positions_d_  = this->positions_ ;
    this->velocities_d_ = this->velocities_;
    this->forces_d_     = this->forces_    ;

    if(this->velocity_initialized_)
    {
        MJOLNIR_LOG_NOTICE(
            "velocity is already given, nothing to initialize in System");
        return ;
    }
    assert(this->has_attribute("temperature"));

    detail::cuda_system_initialize(*this, rng);
    return;
}

template<>
void System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::initialize(rng_type& rng)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // copy initial conditions from CPU to GPU
    this->masses_d_     = this->masses_    ;
    this->rmasses_d_    = this->rmasses_   ;
    this->positions_d_  = this->positions_ ;
    this->velocities_d_ = this->velocities_;
    this->forces_d_     = this->forces_    ;

    if(this->velocity_initialized_)
    {
        MJOLNIR_LOG_NOTICE(
            "velocity is already given, nothing to initialize in System");
        return ;
    }
    assert(this->has_attribute("temperature"));

    detail::cuda_system_initialize(*this, rng);
    return;
}

template<>
void System<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::initialize(rng_type& rng)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // copy initial conditions from CPU to GPU
    this->masses_d_     = this->masses_    ;
    this->rmasses_d_    = this->rmasses_   ;
    this->positions_d_  = this->positions_ ;
    this->velocities_d_ = this->velocities_;
    this->forces_d_     = this->forces_    ;

    if(this->velocity_initialized_)
    {
        MJOLNIR_LOG_NOTICE(
            "velocity is already given, nothing to initialize in System");
        return ;
    }
    assert(this->has_attribute("temperature"));

    detail::cuda_system_initialize(*this, rng);
    return;
}

template<>
void System<CUDASimulatorTraits<double, UnlimitedBoundary>>::sync_configurations()
{
    // assuming that masses will not be changed ...
    this->positions_  = this->positions_d_ ;
    this->velocities_ = this->velocities_d_;
    this->forces_     = this->forces_d_    ;
    return;
}

template<>
void System<CUDASimulatorTraits<float, UnlimitedBoundary>>::sync_configurations()
{
    // assuming that masses will not be changed ...
    this->positions_  = this->positions_d_ ;
    this->velocities_ = this->velocities_d_;
    this->forces_     = this->forces_d_    ;
    return;
}

template<>
void System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::sync_configurations()
{
    // assuming that masses will not be changed ...
    this->positions_  = this->positions_d_ ;
    this->velocities_ = this->velocities_d_;
    this->forces_     = this->forces_d_    ;
    return;
}

template<>
void System<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::sync_configurations()
{
    // assuming that masses will not be changed ...
    this->positions_  = this->positions_d_ ;
    this->velocities_ = this->velocities_d_;
    this->forces_     = this->forces_d_    ;
    return;
}

// ----------------------------------------------------------------------------
// any constructors and `operator=`s are the same as the one that will be
// generated by `= default`, but it should be here because those uses
// `thrust::device_vector`s copy/move constructors and `operator=`s that uses
// device instructions.
//     Since all the implementations are the same, it uses macro.

#define GENERATE_SYSTEM_CONSTRUCTOR_INSTANCES(REAL, BOUNDARY)                   \
    template<>                                                                  \
    System<CUDASimulatorTraits<REAL, BOUNDARY>>::System(                        \
            const std::size_t num_particles, const boundary_type& bd)           \
        : velocity_initialized_(false), boundary_(bd), topology_(num_particles),\
          attributes_(), num_particles_(num_particles),                         \
          masses_     (num_particles), rmasses_     (num_particles),            \
          positions_  (num_particles), velocities_  (num_particles),            \
          forces_     (num_particles),                                          \
          masses_d_   (num_particles), rmasses_d_   (num_particles),            \
          positions_d_(num_particles), velocities_d_(num_particles),            \
          forces_d_   (num_particles)                                           \
    {} /**/

GENERATE_SYSTEM_CONSTRUCTOR_INSTANCES(double, UnlimitedBoundary)
GENERATE_SYSTEM_CONSTRUCTOR_INSTANCES(float,  UnlimitedBoundary)
GENERATE_SYSTEM_CONSTRUCTOR_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_SYSTEM_CONSTRUCTOR_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_SYSTEM_CONSTRUCTOR_INSTANCES

template<> System<CUDASimulatorTraits<double, UnlimitedBoundary>>::~System(){}
template<> System<CUDASimulatorTraits<float,  UnlimitedBoundary>>::~System(){}
template<> System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::~System(){}
template<> System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::~System(){}

#define GENERATE_SYSTEM_COPY_CONSTRUCTOR_INSTANCES(REAL, BOUNDARY)              \
    template<>                                                                  \
    System<CUDASimulatorTraits<REAL, BOUNDARY>>::System(const System& other)    \
        : velocity_initialized_(other.velocity_initialized_),                   \
          boundary_(other.boundary_), topology_(other.topology_),               \
          attributes_(other.attributes_), num_particles_(other.num_particles_), \
          masses_     (other.masses_     ), rmasses_     (other.rmasses_     ), \
          positions_  (other.positions_  ), velocities_  (other.velocities_  ), \
          forces_     (other.forces_     ),                                     \
          masses_d_   (other.masses_d_   ), rmasses_d_   (other.rmasses_d_   ), \
          positions_d_(other.positions_d_), velocities_d_(other.velocities_d_), \
          forces_d_   (other.forces_d_   )                                      \
    {} /**/

GENERATE_SYSTEM_COPY_CONSTRUCTOR_INSTANCES(double, UnlimitedBoundary)
GENERATE_SYSTEM_COPY_CONSTRUCTOR_INSTANCES(float,  UnlimitedBoundary)
GENERATE_SYSTEM_COPY_CONSTRUCTOR_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_SYSTEM_COPY_CONSTRUCTOR_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_SYSTEM_COPY_CONSTRUCTOR_INSTANCES

#define GENERATE_SYSTEM_MOVE_CONSTRUCTOR_INSTANCES(REAL, BOUNDARY)                                     \
    template<>                                                                                         \
    System<CUDASimulatorTraits<REAL, BOUNDARY>>::System(System&& other)                                \
        : velocity_initialized_(other.velocity_initialized_),                                          \
          boundary_   (std::move(other.boundary_)),    topology_     (std::move(other.topology_)),     \
          attributes_ (std::move(other.attributes_ )), num_particles_(std::move(other.num_particles_)),\
          masses_     (std::move(other.masses_     )), rmasses_      (std::move(other.rmasses_      )),\
          positions_  (std::move(other.positions_  )), velocities_   (std::move(other.velocities_   )),\
          forces_     (std::move(other.forces_     )),                                                 \
          masses_d_   (std::move(other.masses_d_   )), rmasses_d_    (std::move(other.rmasses_d_    )),\
          positions_d_(std::move(other.positions_d_)), velocities_d_ (std::move(other.velocities_d_ )),\
          forces_d_   (std::move(other.forces_d_   ))                                                  \
    {} /**/

GENERATE_SYSTEM_MOVE_CONSTRUCTOR_INSTANCES(double, UnlimitedBoundary)
GENERATE_SYSTEM_MOVE_CONSTRUCTOR_INSTANCES(float,  UnlimitedBoundary)
GENERATE_SYSTEM_MOVE_CONSTRUCTOR_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_SYSTEM_MOVE_CONSTRUCTOR_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_SYSTEM_MOVE_CONSTRUCTOR_INSTANCES

#define GENERATE_SYSTEM_COPY_ASSIGNER_INSTANCES(REAL, BOUNDARY)                 \
    template<>                                                                  \
    System<CUDASimulatorTraits<REAL, BOUNDARY>>&                                \
    System<CUDASimulatorTraits<REAL, BOUNDARY>>::operator=(const System& other) \
    {                                                                           \
        this->velocity_initialized_ = other.velocity_initialized_;              \
        this->boundary_      = other.boundary_     ;                            \
        this->topology_      = other.topology_     ;                            \
        this->attributes_    = other.attributes_   ;                            \
        this->num_particles_ = other.num_particles_;                            \
        this->masses_        = other.masses_       ;                            \
        this->rmasses_       = other.rmasses_      ;                            \
        this->positions_     = other.positions_    ;                            \
        this->velocities_    = other.velocities_   ;                            \
        this->forces_        = other.forces_       ;                            \
        this->masses_d_      = other.masses_d_     ;                            \
        this->rmasses_d_     = other.rmasses_d_    ;                            \
        this->positions_d_   = other.positions_d_  ;                            \
        this->velocities_d_  = other.velocities_d_ ;                            \
        this->forces_d_      = other.forces_d_     ;                            \
        return *this;                                                           \
    } /**/

GENERATE_SYSTEM_COPY_ASSIGNER_INSTANCES(double, UnlimitedBoundary)
GENERATE_SYSTEM_COPY_ASSIGNER_INSTANCES(float,  UnlimitedBoundary)
GENERATE_SYSTEM_COPY_ASSIGNER_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_SYSTEM_COPY_ASSIGNER_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_SYSTEM_COPY_ASSIGNER_INSTANCES

#define GENERATE_SYSTEM_MOVE_ASSIGNER_INSTANCES(REAL, BOUNDARY)            \
    template<>                                                             \
    System<CUDASimulatorTraits<REAL, BOUNDARY>>&                           \
    System<CUDASimulatorTraits<REAL, BOUNDARY>>::operator=(System&& other) \
    {                                                                      \
        this->velocity_initialized_ = other.velocity_initialized_;         \
        this->boundary_      = std::move(other.boundary_     );            \
        this->topology_      = std::move(other.topology_     );            \
        this->attributes_    = std::move(other.attributes_   );            \
        this->num_particles_ = std::move(other.num_particles_);            \
        this->masses_        = std::move(other.masses_       );            \
        this->rmasses_       = std::move(other.rmasses_      );            \
        this->positions_     = std::move(other.positions_    );            \
        this->velocities_    = std::move(other.velocities_   );            \
        this->forces_        = std::move(other.forces_       );            \
        this->masses_d_      = std::move(other.masses_d_     );            \
        this->rmasses_d_     = std::move(other.rmasses_d_    );            \
        this->positions_d_   = std::move(other.positions_d_  );            \
        this->velocities_d_  = std::move(other.velocities_d_ );            \
        this->forces_d_      = std::move(other.forces_d_     );            \
        return *this;                                                      \
    } /**/

GENERATE_SYSTEM_MOVE_ASSIGNER_INSTANCES(double, UnlimitedBoundary)
GENERATE_SYSTEM_MOVE_ASSIGNER_INSTANCES(float,  UnlimitedBoundary)
GENERATE_SYSTEM_MOVE_ASSIGNER_INSTANCES(double, CuboidalPeriodicBoundary)
GENERATE_SYSTEM_MOVE_ASSIGNER_INSTANCES(float,  CuboidalPeriodicBoundary)

#undef GENERATE_SYSTEM_MOVE_ASSIGNER_INSTANCES

// ---------------------------------------------------------------------------
// class template instanciation

template class System<CUDASimulatorTraits<double, UnlimitedBoundary>>;
template class System<CUDASimulatorTraits<float,  UnlimitedBoundary>>;
template class System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
