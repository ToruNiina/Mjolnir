#ifndef MJOLNIR_CUDA_BAOAB_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CUDA_BAOAB_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/cuda/CUDASimulatorTraits.hpp>
#include <mjolnir/cuda/System.hpp>
#include <mjolnir/cuda/RandomNumberGenerator.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace mjolnir
{

// a specialization of BAOAB Langevin integrator for OpenMP implementation
template<typename realT, template<typename, typename> class boundaryT>
class BAOABLangevinIntegrator<CUDASimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = CUDASimulatorTraits<realT, boundaryT>;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = ForceField<traits_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;

  public:

    BAOABLangevinIntegrator(const real_type dt,
                            const std::vector<real_type>& gamma);

    BAOABLangevinIntegrator(const BAOABLangevinIntegrator&);
    BAOABLangevinIntegrator(BAOABLangevinIntegrator&&);
    BAOABLangevinIntegrator& operator=(const BAOABLangevinIntegrator&);
    BAOABLangevinIntegrator& operator=(BAOABLangevinIntegrator&&);
    ~BAOABLangevinIntegrator();

    void initialize(system_type& sys, forcefield_type& ff, rng_type&);

    real_type step(const real_type  time, system_type& sys,
                   forcefield_type& ff,   rng_type&    rng);

    void update(const system_type& sys)
    {
        if(!sys.has_attribute("temperature"))
        {
            throw std::out_of_range("mjolnir::BAOABLangevinIntegrator: "
                "Langevin Integrator requires reference temperature, but "
                "`temperature` is not found in `system.attribute`.");
        }
        this->temperature_ = sys.attribute("temperature");
        this->reset_parameters(sys);
        return;
    }

    real_type delta_t()     const noexcept {return dt_;}
    real_type halfdt()      const noexcept {return halfdt_;}
    real_type temperature() const noexcept {return temperature_;}

    thrust::host_vector<real_type> const& parameters() const noexcept
    {
        return gammas_h_;
    }

    thrust::device_vector<real_type>& gammas_device()        noexcept {return this->gammas_;}
    thrust::device_vector<real_type>& exp_gamma_dt_device()  noexcept {return this->exp_gamma_dt_;}
    thrust::device_vector<real_type>& noise_coeff_device()   noexcept {return this->noise_coeff_; }
    thrust::device_vector<real_type>& random_forces_device() noexcept {return this->random_forces_;}

    thrust::device_vector<real_type> const& gammas_device()        const noexcept {return this->gammas_;}
    thrust::device_vector<real_type> const& exp_gamma_dt_device()  const noexcept {return this->exp_gamma_dt_;}
    thrust::device_vector<real_type> const& noise_coeff_device()   const noexcept {return this->noise_coeff_; }
    thrust::device_vector<real_type> const& random_forces_device() const noexcept {return this->random_forces_;}

  private:

    void reset_parameters(const system_type& sys);

  private:
    real_type dt_;
    real_type halfdt_;
    real_type temperature_;

    thrust::host_vector<real_type>   gammas_h_;
    thrust::device_vector<real_type> gammas_;
    thrust::device_vector<real_type> exp_gamma_dt_;
    thrust::device_vector<real_type> noise_coeff_;
    thrust::device_vector<real_type> random_forces_;
};

// ---------------------------------------------------------------------------
// declaration of BAOABLangevinIntegrator<CUDA>::step();

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng);

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng);

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng);

template<>
typename BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::real_type
BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::step(
    const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng);

// ---------------------------------------------------------------------------
// declaration of BAOABLangevinIntegrator<CUDA>::initialize();

template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&);
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&);
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&);
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>>::initialize(
        system_type& sys, forcefield_type& ff, rng_type&);

// ---------------------------------------------------------------------------
// definitions of BAOABLangevinIntegrator<CUDA>::reset_parameters();

template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>
    >::reset_parameters(const system_type& sys);
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, UnlimitedBoundary>
    >::reset_parameters(const system_type& sys);
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>
    >::reset_parameters(const system_type& sys);
template<>
void BAOABLangevinIntegrator<CUDASimulatorTraits<float, CuboidalPeriodicBoundary>
    >::reset_parameters(const system_type& sys);

// ---------------------------------------------------------------------------
// declaration of BAOABLangevinIntegrator<CUDA>::constructors and others.
//
// These uses device_vector, so use device functions. If it would be generated
// through .cpp file, it causes compilation error. Thus those should be defined
// in .cu files, and it requires declaration to prohibit automatic generation.

template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >::BAOABLangevinIntegrator(const real_type, const std::vector<real_type>&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::BAOABLangevinIntegrator(const real_type, const std::vector<real_type>&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::BAOABLangevinIntegrator(const real_type, const std::vector<real_type>&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::BAOABLangevinIntegrator(const real_type, const std::vector<real_type>&);

template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >::~BAOABLangevinIntegrator();
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::~BAOABLangevinIntegrator();
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::~BAOABLangevinIntegrator();
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::~BAOABLangevinIntegrator();

template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >::BAOABLangevinIntegrator(const BAOABLangevinIntegrator&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::BAOABLangevinIntegrator(const BAOABLangevinIntegrator&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::BAOABLangevinIntegrator(const BAOABLangevinIntegrator&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::BAOABLangevinIntegrator(const BAOABLangevinIntegrator&);

template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >::BAOABLangevinIntegrator(BAOABLangevinIntegrator&&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::BAOABLangevinIntegrator(BAOABLangevinIntegrator&&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::BAOABLangevinIntegrator(BAOABLangevinIntegrator&&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::BAOABLangevinIntegrator(BAOABLangevinIntegrator&&);

template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >& BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >::operator=(const BAOABLangevinIntegrator&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>& BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::operator=(const BAOABLangevinIntegrator&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >& BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::operator=(const BAOABLangevinIntegrator&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>& BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::operator=(const BAOABLangevinIntegrator&);

template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >& BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>       >::operator=(BAOABLangevinIntegrator&&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>& BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>::operator=(BAOABLangevinIntegrator&&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >& BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>       >::operator=(BAOABLangevinIntegrator&&);
template<> BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>& BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>::operator=(BAOABLangevinIntegrator&&);

extern template class BAOABLangevinIntegrator<CUDASimulatorTraits<double, UnlimitedBoundary>>;
extern template class BAOABLangevinIntegrator<CUDASimulatorTraits<float,  UnlimitedBoundary>>;
extern template class BAOABLangevinIntegrator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class BAOABLangevinIntegrator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif /* MJOLNIR_BAOAB_LANGEVIN_INTEGRATOR */
