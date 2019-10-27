set(mjolnir_omp_nondet_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/RandomNumberGenerator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BAOABLangevinIntegrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/UnderdampedLangevinIntegrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairLennardJonesInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/UnlimitedGridCellList.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGridCellList.cpp"
    )
set(mjolnir_cpp_files
    ${mjolnir_cpp_files}
    ${mjolnir_omp_nondet_cpp_files} PARENT_SCOPE)
