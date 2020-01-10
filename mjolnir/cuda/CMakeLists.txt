# compile cuda stuffs as a library

set(mjolnir_cuda_cu_files
    "${CMAKE_CURRENT_SOURCE_DIR}/BAOABLangevinIntegrator.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/RandomNumberGenerator.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/System.cu"
    )

add_library(mjolnir_cuda STATIC ${mjolnir_cuda_cu_files})
set_target_properties(mjolnir_cuda PROPERTIES
    COMPILE_FLAGS "-ccbin ${CMAKE_CXX_COMPILER} -std=c++11 -O3 --use_fast_math --expt-extended-lambda --expt-relaxed-constexpr")

# add cuda specific cpp files

set(mjolnir_cuda_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/ObserverContainer.cpp"
    )

set(mjolnir_cpp_files
    ${mjolnir_cpp_files}
    ${mjolnir_cuda_cpp_files} PARENT_SCOPE)
