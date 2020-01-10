#ifndef MJOLNIR_UTIL_MACRO_HPP
#define MJOLNIR_UTIL_MACRO_HPP

// MJOLNIR_STRINGIZE(hoge) -> "hoge"
#define MJOLNIR_STRINGIZE_AUX(x) #x
#define MJOLNIR_STRINGIZE(x)     MJOLNIR_STRINGIZE_AUX(x)

// MJOLNIR_FUNC_NAME, expanded into function name.
#if defined(__GNUC__)
#  define MJOLNIR_FUNC_NAME __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
#  define MJOLNIR_FUNC_NAME __FUNCSIG__
#else
#  define MJOLNIR_FUNC_NAME __func__
#endif

// for CUDA. if the compiler is not nvcc, these macros will be expanded to empty
#ifndef MJOLNIR_CUDA_HOST_DEVICE
#  ifdef __NVCC__
#    define MJOLNIR_CUDA_HOST_ONLY   __host__
#    define MJOLNIR_CUDA_DEVICE_ONLY __device__
#    define MJOLNIR_CUDA_HOST_DEVICE __host__ __device__
#    define MJOLNIR_CUDA_FORCEINLINE __forceinline__
#  else
#    define MJOLNIR_CUDA_HOST_ONLY
#    define MJOLNIR_CUDA_DEVICE_ONLY
#    define MJOLNIR_CUDA_HOST_DEVICE
#    define MJOLNIR_CUDA_FORCEINLINE
#  endif
#endif// MJOLNIR_CUDA_HOST_DEVICE

#endif // MJOLNIR_UTIL_MACRO_HPP
