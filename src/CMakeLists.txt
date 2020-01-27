# ----------------------------------------------------------------------------
# add executables

# message(STATUS ${mjolnir_source_files})

add_executable(mjolnir mjolnir.cpp)
set_target_properties(mjolnir PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    COMPILE_FLAGS "${MJOLNIR_WARNING_FLAGS} ${MJOLNIR_OPTIMIZATON_FLAGS}")

if(SEPARATE_BUILD)
    add_library(mjolnir_core STATIC ${mjolnir_source_files})
    set_target_properties(mjolnir_core PROPERTIES
        COMPILE_FLAGS "${MJOLNIR_WARNING_FLAGS} ${MJOLNIR_OPTIMIZATON_FLAGS}")
    target_link_libraries(mjolnir mjolnir_core)

    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0")
        target_precompile_headers(mjolnir_core
            PRIVATE "${PROJECT_SOURCE_DIR}/extlib/toml/toml.hpp")

        target_precompile_headers(mjolnir REUSE_FROM mjolnir_core)
    endif()
endif()

if(OpenMP_FOUND AND USE_OPENMP)
    message(STATUS "adding OpenMP flags ${OpenMP_CXX_FLAGS} to mjolnir ...")
    target_link_libraries(mjolnir ${OpenMP_CXX_LIBRARIES})

    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        # After CMake 3.13, we can use target_link_options
        set_target_properties(mjolnir PROPERTIES LINK_FLAGS "-parallel")
    endif()
else()
    message(STATUS "Ignoring OpenMP ...")
endif()
