#ifndef MJOLNIR_CUDA_CUDA_SIMULATOR_TRAITS
#define MJOLNIR_CUDA_CUDA_SIMULATOR_TRAITS
#include <mjolnir/cuda/Vector.hpp>

namespace mjolnir
{

// full CUDA

template<typename realT, template<typename, typename> class boundaryT>
struct CUDASimulatorTraits
{
    using real_type       = realT;
    using coordinate_type = typename cuda::vector_type_of<real_type>::type;

//     using matrix33_type = math::Matrix<real_type, 3, 3>;
//     using matrix44_type = math::Matrix<real_type, 4, 4>;

    using boundary_type = boundaryT<real_type, coordinate_type>;
};

template<typename T>
struct is_cuda_simulator_traits : std::false_type{};

template<typename realT, template<typename, typename> class boundaryT>
struct is_cuda_simulator_traits<CUDASimulatorTraits<realT, boundaryT>>: std::true_type{};


} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
