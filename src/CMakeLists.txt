add_library(demsilib
    demsi_partition.cpp
    demsi_logging.cpp
    demsi_kokkosparser.cpp
    demsi_grid.cpp
    demsi_grid_io.cpp
    demsi_lmp_instance.cpp
    demsi_configs.cpp
    demsi_tinyxml2.cpp
    demsi_time.cpp
    demsi_particles.cpp
    demsi_particles_io.cpp
    demsi_forcing.cpp
    demsi_external_force.cpp
    demsi_file_utils.cpp
    demsi_interpolation.cpp
    demsi_contacts.cpp
    demsi_column.cpp
    demsi_column_icepack.F90
    demsi_column_tracers.cpp
    demsi_column_variables.cpp
    demsi_column_transfer.cpp
    demsi_communication.cpp
    demsi_polygon.cpp
    demsi_initial_tessellation.cpp
    demsi_initial_tessellation_io.cpp
    demsi_remapping.cpp
    demsi_remapping_flux.cpp
    demsi_remapping_obr.cpp
    demsi_remapping_debug.cpp
    demsi_ocean.cpp
    demsi_ridging.cpp
    demsi_ridging.F90
    demsi_diagnostics.cpp
    )

# RPATH must be appended before the executable requiring that entry from RPATH is declared
set(CMAKE_INSTALL_RPATH ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_INSTALL_RPATH})
add_executable(demsi
    demsi.cpp
    )

#make it show up in the root directory
set_target_properties(demsi PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

message("-- CMAKE_SHARED_LIBRARY_SUFFIX: ${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(LAMMPS_INCLUDE_DIR "${LAMMPS_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
message("-- LAMMPS include: ${LAMMPS_INCLUDE_DIR}")
set(LAMMPS_LIB "${LAMMPS_PREFIX}/${CMAKE_INSTALL_LIBDIR}/liblammps${CMAKE_SHARED_LIBRARY_SUFFIX}")
message("-- LAMMPS lib: ${LAMMPS_LIB}")
set(KOKKOS_INCLUDE_DIR "${LAMMPS_PREFIX}/include/kokkos")
message("-- KOKKOS include: ${KOKKOS_INCLUDE_DIR}")
set(KOKKOS_LIB "${LAMMPS_PREFIX}/lib/libkokkos${CMAKE_SHARED_LIBRARY_SUFFIX}")
message("-- KOKKOS lib: ${KOKKOS_LIB}")
set(petsc_INCLUDE_DIR "${petsc_PREFIX}/include")
message("-- petsc include: ${petsc_INCLUDE_DIR}")
set(petsc_LIB "${petsc_PREFIX}/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}")
message("-- petsc lib: ${petsc_LIB}")
set(permon_INCLUDE_DIR "${permon_PREFIX}/include")
message("-- permon include: ${permon_INCLUDE_DIR}")
set(permon_LIB "${permon_PREFIX}/lib/libpermon${CMAKE_SHARED_LIBRARY_SUFFIX}")
message("-- permon lib: ${permon_LIB}")

add_dependencies(demsilib LAMMPS)
target_include_directories(demsilib PUBLIC "${LAMMPS_INCLUDE_DIR}")
target_link_libraries(demsilib PUBLIC "${LAMMPS_LIB}")
target_include_directories(demsilib PUBLIC "${KOKKOS_INCLUDE_DIR}")
target_link_libraries(demsilib PUBLIC "${KOKKOS_LIB}")
target_link_libraries(demsilib PUBLIC Icepack)
target_include_directories(demsilib PUBLIC "${netCDF_INCLUDE_DIR}")
target_link_libraries(demsilib PUBLIC "${netCDF_LIBRARY}")
target_include_directories(demsilib PUBLIC "${petsc_INCLUDE_DIR}")
target_link_libraries(demsilib PUBLIC "${petsc_LIB}")
target_include_directories(demsilib PUBLIC "${permon_INCLUDE_DIR}")
target_link_libraries(demsilib PUBLIC "${permon_LIB}")

target_link_libraries(demsi PUBLIC demsilib)
