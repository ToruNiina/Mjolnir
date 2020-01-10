#ifndef MJOLNIR_CUDA_SYSTEM_HPP
#define MJOLNIR_CUDA_SYSTEM_HPP
#include <mjolnir/cuda/CUDASimulatorTraits.hpp>
#include <mjolnir/cuda/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace mjolnir
{

// System<CUDA> manages memory on both CPU and GPU to store positions,
// velocities, forces and masses.
// In this CUDA implementation, everything will be calculated on GPU and
// only when it outputs current state to files, it synchronizes. Therefore,
// unlike System<OpenMP>, System<CUDA> does not have any member method to
// sync and merge forces calculated on CPU and GPU. So TAKE CARE when you
// access to a value via access methods.
template<typename realT, template<typename, typename> class boundaryT>
class System<CUDASimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = CUDASimulatorTraits<realT, boundaryT>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using boundary_type   = typename traits_type::boundary_type;
    using topology_type   = Topology;
    using attribute_type  = std::map<std::string, real_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;

    using string_type              = std::string;
    using particle_type            = Particle<real_type, coordinate_type>;
    using particle_view_type       = ParticleView<real_type, coordinate_type>;
    using particle_const_view_type = ParticleConstView<real_type, coordinate_type>;

  public:

    // Since most of the function uses CUDA device stuff, they should be defined
    // in `.cu` file. Otherwise, compilation fails.

    System(const std::size_t num_particles, const boundary_type& bd);
    System(const System&);
    System(System&&);
    System& operator=(const System&);
    System& operator=(System&&);
    ~System();

    void initialize(rng_type& rng);

    coordinate_type adjust_direction(coordinate_type dr) const noexcept
    {return boundary_.adjust_direction(dr);}
    coordinate_type  adjust_position(coordinate_type dr) const noexcept
    {return boundary_.adjust_position(dr);}

    std::size_t size() const noexcept {return num_particles_;}

    particle_view_type operator[](std::size_t i) noexcept
    {
        return particle_view_type{
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_[i],
            this->topology_.name_of (i, std::nothrow),
            this->topology_.group_of(i, std::nothrow)
        };
    }
    particle_const_view_type operator[](std::size_t i) const noexcept
    {
        return particle_const_view_type{
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_[i],
            this->topology_.name_of (i, std::nothrow),
            this->topology_.group_of(i, std::nothrow)
        };
    }

    template<typename Vec>
    void check_size(const Vec& v, const std::size_t i, const char* name) const
    {
        if(v.size() <= i)
        {
            throw_exception<std::out_of_range>("System::at(", i, "), ", name,
                    " has size = ", v.size());
        }
    }

    particle_view_type at(std::size_t i)
    {
        // thrust::host_vector does not have at member method
        this->check_size(masses_    , i, "masses");
        this->check_size(rmasses_   , i, "rmasses");
        this->check_size(positions_ , i, "positions");
        this->check_size(velocities_, i, "velocities");
        this->check_size(forces_    , i, "forces");

        return particle_view_type{
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_[i],
            this->topology_.name_of (i),
            this->topology_.group_of(i)
        };
    }
    particle_const_view_type at(std::size_t i) const
    {
        this->check_size(masses_    , i, "masses");
        this->check_size(rmasses_   , i, "rmasses");
        this->check_size(positions_ , i, "positions");
        this->check_size(velocities_, i, "velocities");
        this->check_size(forces_    , i, "forces");

        return particle_const_view_type{
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_[i],
            this->topology_.name_of (i),
            this->topology_.group_of(i)
        };
    }

    real_type  mass (std::size_t i) const noexcept {return masses_[i];}
    real_type& mass (std::size_t i)       noexcept {return masses_[i];}
    real_type  rmass(std::size_t i) const noexcept {return rmasses_[i];}
    real_type& rmass(std::size_t i)       noexcept {return rmasses_[i];}

    coordinate_type const& position(std::size_t i) const noexcept {return positions_[i];}
    coordinate_type&       position(std::size_t i)       noexcept {return positions_[i];}
    coordinate_type const& velocity(std::size_t i) const noexcept {return velocities_[i];}
    coordinate_type&       velocity(std::size_t i)       noexcept {return velocities_[i];}
    coordinate_type const& force   (std::size_t i) const noexcept {return forces_[i];}
    coordinate_type&       force   (std::size_t i)       noexcept {return forces_[i];}

    // -----------------------------------------------------------------------
    // host data access

    // As mentioned at the top of this file, this implementation does not use
    // CPU. So this function always overwrite CPU value by copying GPU values.
    void sync_configurations();

    thrust::host_vector<real_type> const& masses_host()  const noexcept {return masses_;}
    thrust::host_vector<real_type> &      masses_host()        noexcept {return masses_;}
    thrust::host_vector<real_type> const& rmasses_host() const noexcept {return rmasses_;}
    thrust::host_vector<real_type> &      rmasses_host()       noexcept {return rmasses_;}

    thrust::host_vector<coordinate_type> const& positions_host()  const noexcept {return positions_;}
    thrust::host_vector<coordinate_type> &      positions_host()        noexcept {return positions_;}
    thrust::host_vector<coordinate_type> const& velocities_host() const noexcept {return velocities_;}
    thrust::host_vector<coordinate_type> &      velocities_host()       noexcept {return velocities_;}
    thrust::host_vector<coordinate_type> const& forces_host()     const noexcept {return forces_;}
    thrust::host_vector<coordinate_type> &      forces_host()           noexcept {return forces_;}

    // -----------------------------------------------------------------------
    // device data access

    thrust::device_vector<real_type> const& masses_device()  const noexcept {return masses_d_;}
    thrust::device_vector<real_type> &      masses_device()        noexcept {return masses_d_;}
    thrust::device_vector<real_type> const& rmasses_device() const noexcept {return rmasses_d_;}
    thrust::device_vector<real_type> &      rmasses_device()       noexcept {return rmasses_d_;}

    thrust::device_vector<coordinate_type> const& positions_device()  const noexcept {return positions_d_;}
    thrust::device_vector<coordinate_type> &      positions_device()        noexcept {return positions_d_;}
    thrust::device_vector<coordinate_type> const& velocities_device() const noexcept {return velocities_d_;}
    thrust::device_vector<coordinate_type> &      velocities_device()       noexcept {return velocities_d_;}
    thrust::device_vector<coordinate_type> const& forces_device()     const noexcept {return forces_d_;}
    thrust::device_vector<coordinate_type> &      forces_device()           noexcept {return forces_d_;}

    // -----------------------------------------------------------------------
    // other stuff

    string_type const& name (std::size_t i) const noexcept {return topology_.name_of (i, std::nothrow);}
    string_type&       name (std::size_t i)       noexcept {return topology_.name_of (i, std::nothrow);}
    string_type const& group(std::size_t i) const noexcept {return topology_.group_of(i, std::nothrow);}
    string_type&       group(std::size_t i)       noexcept {return topology_.group_of(i, std::nothrow);}

    boundary_type&       boundary()       noexcept {return boundary_;}
    boundary_type const& boundary() const noexcept {return boundary_;}
    topology_type&       topology()       noexcept {return topology_;}
    topology_type const& topology() const noexcept {return topology_;}

    real_type  attribute(const std::string& key) const {return attributes_.at(key);}
    real_type& attribute(const std::string& key)       {return attributes_[key];}
    bool   has_attribute(const std::string& key) const {return attributes_.count(key) == 1;}
    attribute_type const& attributes() const noexcept {return attributes_;}

    bool  velocity_initialized() const noexcept {return velocity_initialized_;}
    bool& velocity_initialized()       noexcept {return velocity_initialized_;}

  private:

    bool           velocity_initialized_;
    boundary_type  boundary_;
    topology_type  topology_;
    attribute_type attributes_;
    std::size_t    num_particles_;

    thrust::host_vector<real_type>       masses_;
    thrust::host_vector<real_type>       rmasses_; // 1 / mass
    thrust::host_vector<coordinate_type> positions_;
    thrust::host_vector<coordinate_type> velocities_;
    thrust::host_vector<coordinate_type> forces_;

    // containers on GPU memory
    thrust::device_vector<real_type>       masses_d_;
    thrust::device_vector<real_type>       rmasses_d_;
    thrust::device_vector<coordinate_type> positions_d_;
    thrust::device_vector<coordinate_type> velocities_d_;
    thrust::device_vector<coordinate_type> forces_d_;
};

// ---------------------------------------------------------------------------
// declaration of member functions that use device codes.

template<> void System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::initialize(rng_type& rng);
template<> void System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::initialize(rng_type& rng);
template<> void System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::initialize(rng_type& rng);
template<> void System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::initialize(rng_type& rng);

template<> void System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::sync_configurations();
template<> void System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::sync_configurations();
template<> void System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::sync_configurations();
template<> void System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::sync_configurations();

// ---------------------------------------------------------------------------
// declaration of System<CUDA>::constructors and others.
//
// These uses device_vector, so use device functions. If it would be generated
// through .cpp file, it causes compilation error. Thus those should be defined
// in .cu files, and it requires declaration to prohibit automatic generation.

template<> System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::System(const std::size_t num_particles, const boundary_type& bd);
template<> System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::System(const std::size_t num_particles, const boundary_type& bd);
template<> System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::System(const std::size_t num_particles, const boundary_type& bd);
template<> System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::System(const std::size_t num_particles, const boundary_type& bd);

template<> System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::System(const System&);
template<> System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::System(const System&);
template<> System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::System(const System&);
template<> System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::System(const System&);

template<> System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::System(System&&);
template<> System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::System(System&&);
template<> System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::System(System&&);
template<> System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::System(System&&);

template<> System<CUDASimulatorTraits<double, UnlimitedBoundary>       >& System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::operator=(const System&);
template<> System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >& System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::operator=(const System&);
template<> System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>& System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::operator=(const System&);
template<> System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>& System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::operator=(const System&);

template<> System<CUDASimulatorTraits<double, UnlimitedBoundary>       >& System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::operator=(System&&);
template<> System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >& System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::operator=(System&&);
template<> System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>& System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::operator=(System&&);
template<> System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>& System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::operator=(System&&);

template<> System<CUDASimulatorTraits<double, UnlimitedBoundary>       >::~System();
template<> System<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::~System();
template<> System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::~System();
template<> System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::~System();

extern template class System<CUDASimulatorTraits<double, UnlimitedBoundary>>;
extern template class System<CUDASimulatorTraits<float,  UnlimitedBoundary>>;
extern template class System<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class System<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif// MJOLNIR_CUDA_SYSTEM_HPP
