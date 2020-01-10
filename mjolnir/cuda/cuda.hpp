#ifndef MJOLNIR_CUDA_CUDA_HPP
#define MJOLNIR_CUDA_CUDA_HPP

// This file is a meta-header file that just includes everything that are needed
// to use CUDA implementation. This file is introduced in order to make it
// simple to turn on/off CUDA support by preprocessor macro, like in the
// following way.
//
// ```cpp
// #ifdef WITH_CUDA
// #include <mjolnir/cuda/cuda.hpp>
// #endif
// ```

#ifndef MJOLNIR_SEPARATE_BUILD
#error "cuda support requires SEPARATE_BUILD=ON"
#endif // MJOLNIR_SEPARATE_BUILD

#include <mjolnir/cuda/CUDASimulatorTraits.hpp>
#include <mjolnir/cuda/System.hpp>
#include <mjolnir/cuda/RandomNumberGenerator.hpp>
#include <mjolnir/cuda/BAOABLangevinIntegrator.hpp>
#include <mjolnir/cuda/ObserverContainer.hpp>

#endif// MJOLNIR_CUDA_CUDA_HPP
