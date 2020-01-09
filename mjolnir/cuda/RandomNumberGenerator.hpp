#ifndef MJOLNIR_CUDA_RANDOM_NUMBER_GENERATOR_HPP
#define MJOLNIR_CUDA_RANDOM_NUMBER_GENERATOR_HPP
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/cuda/CUDASimulatorTraits.hpp>
#include <thrust/device_vector.h>
#include <cuda.h>
#include <curand.h>
#include <random>

namespace mjolnir
{
// wrapper functions to call appropreate cuRAND API by overloading double/float
namespace detail
{
curandStatus_t curand_generate_uniform_01(curandGenerator_t gen, thrust::device_vector<double>& buf);
curandStatus_t curand_generate_uniform_01(curandGenerator_t gen, thrust::device_vector<float >& buf);
curandStatus_t curand_generate_gaussian  (curandGenerator_t gen, thrust::device_vector<double>& buf);
curandStatus_t curand_generate_gaussian  (curandGenerator_t gen, thrust::device_vector<float >& buf);
} // detail

template<typename realT, template<typename, typename> class boundaryT>
class RandomNumberGenerator<CUDASimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = CUDASimulatorTraits<realT, boundaryT>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

  public:

    explicit RandomNumberGenerator(const std::uint32_t seed)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // initialize cuRAND generator
        {
            const curandStatus_t stat = curandCreateGenerator(
                    std::addressof(this->generator_),
                    CURAND_RNG_PSEUDO_DEFAULT); // XORWOW is the current default
            if(stat != CURAND_STATUS_SUCCESS)
            {
                MJOLNIR_LOG_ERROR("curandCreateGenerator failed with status ", stat);
                curandDestroyGenerator(this->generator_);
                throw std::runtime_error("curandCreateGenerator failed");
            }
        }
        // set seed.
        //
        // To use rng both on CPU and GPU, generate 2 seeds from the seed passed.
        std::seed_seq sseq{seed};
        std::array<std::uint32_t, 2> buf;
        sseq.generate(buf.begin(), buf.end());

        // set CPU rng.
        this->rng_.seed(buf.front());
        {
            // set GPU rng.
            const curandStatus_t stat = curandSetPseudoRandomGeneratorSeed(
                    this->generator_, buf.back());
            if(stat != CURAND_STATUS_SUCCESS)
            {
                MJOLNIR_LOG_ERROR("curandSetPseudoRandomGeneratorSeed failed "
                                  "with status ", stat);
                curandDestroyGenerator(this->generator_);
                throw std::runtime_error(
                        "curandSetPseudoRandomGeneratorSeed failed");
            }
        }
    }
    ~RandomNumberGenerator() noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        const curandStatus_t stat = curandDestroyGenerator(generator_);
        if(stat != CURAND_STATUS_SUCCESS)
        {
            MJOLNIR_LOG_ERROR("curandDestroyGenerator failed with status ", stat);
        }
    }

    void fill_with_uniform_01(thrust::device_vector<real_type>& buf)
    {
        if(buf.size() == 0) {return ;}

        const auto stat = detail::curand_generate_uniform_01(generator_, buf);
        assert(stat == CURAND_STATUS_SUCCESS);
        (void)stat;
        return;
    }
    void fill_with_gaussian(thrust::device_vector<real_type>& buf)
    {
        if(buf.size() == 0) {return ;}

        const auto stat = detail::curand_generate_gaussian(generator_, buf);
        assert(stat == CURAND_STATUS_SUCCESS);
        (void)stat;
        return;
    }

    real_type uniform_real01()
    {
        return std::generate_canonical<
            real_type, std::numeric_limits<real_type>::digits
            >(this->rng_);
    }
    real_type uniform_real(const real_type min, const real_type max)
    {
        return this->uniform_real01() * (max - min) + min;
    }

    real_type gaussian()
    {
        return this->nrm_(this->rng_);
    }
    real_type gaussian(const real_type mean, const real_type stddev)
    {
        return this->nrm_(this->rng_) * stddev + mean;
    }


  private:
    std::uint32_t     seed_;
    curandGenerator_t generator_;
    std::mt19937      rng_;
    std::normal_distribution<real_type> nrm_;
};

extern template class RandomNumberGenerator<CUDASimulatorTraits<double, UnlimitedBoundary>>;
extern template class RandomNumberGenerator<CUDASimulatorTraits<float,  UnlimitedBoundary>>;
extern template class RandomNumberGenerator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class RandomNumberGenerator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif /*MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR*/
