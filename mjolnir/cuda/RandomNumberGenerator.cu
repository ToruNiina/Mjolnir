#include <mjolnir/cuda/RandomNumberGenerator.hpp>

namespace mjolnir
{
namespace detail
{
curandStatus_t curand_generate_uniform_01(curandGenerator_t gen, thrust::device_vector<double>& buf)
{
    return curandGenerateUniformDouble(gen, buf.data().get(), buf.size());
}
curandStatus_t curand_generate_uniform_01(curandGenerator_t gen, thrust::device_vector<float>& buf)
{
    return curandGenerateUniform(gen, buf.data().get(), buf.size());
}
curandStatus_t curand_generate_gaussian(curandGenerator_t gen, thrust::device_vector<double>& buf)
{
    return curandGenerateNormalDouble(
            gen, buf.data().get(), buf.size(), /*mean*/ 0.0, /*stddev*/ 1.0);
}
curandStatus_t curand_generate_gaussian(curandGenerator_t gen, thrust::device_vector<float>& buf)
{
    return curandGenerateNormal(
            gen, buf.data().get(), buf.size(), /*mean*/ 0.0, /*stddev*/ 1.0);
}
} // detail

// ---------------------------------------------------------------------------
// class template instanciation
//

template class RandomNumberGenerator<CUDASimulatorTraits<double, UnlimitedBoundary>>;
template class RandomNumberGenerator<CUDASimulatorTraits<float,  UnlimitedBoundary>>;
template class RandomNumberGenerator<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class RandomNumberGenerator<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;
}
