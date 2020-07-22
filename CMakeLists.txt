cmake_minimum_required(VERSION ${CMAKE_VERSION})
cmake_policy(VERSION ${CMAKE_VERSION})

project(DEMSI LANGUAGES C CXX Fortran)

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------

#to build LAMMPS along with DEMSI
include(ExternalProject)
#to figure out where LAMMPS is going to put its library
include(GNUInstallDirs)
#to check compiler flags are supported
include(CheckCXXCompilerFlag)

#-------------------------------------------------------------------------------
# C++ compiler flags
#-------------------------------------------------------------------------------

# flags for all compilers and build types
set(CMAKE_CXX_FLAGS "-g --std=c++11")

# check if DEMSI_USE_CUDA set in environment is set to ON
string(TOUPPER "$ENV{DEMSI_USE_CUDA}" USE_CUDA)
if ("${USE_CUDA}" STREQUAL "ON")
  message( STATUS "DEMSI_USE_CUDA set to ON.")
  set(DEMSI_USE_CUDA ON)
else()
  set(DEMSI_USE_CUDA OFF)
endif()

if (DEMSI_USE_OPENMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  message( STATUS "DEMSI_USE_OPENMP set to ON.")
endif()

if (DEMSI_USE_CUDA) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda")
endif()

# release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# debug flags
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

  set(CMAKE_CXX_FLAGS_DEBUG "-O0")

  unset(CXX_COMPILER_SUPPORTS_SANITIZE_BOUNDS CACHE)
  unset(CXX_COMPILER_SUPPORTS_SANITIZE_ADDRESS CACHE)
  unset(CXX_COMPILER_SUPPORTS_SANITIZE_UNDEFINED CACHE)
  unset(CXX_COMPILER_SUPPORTS_SANITIZE_ADDRESS_USE_AFTER_SCOPE CACHE)

  if (NOT DEMSI_USE_CUDA) 
    CHECK_CXX_COMPILER_FLAG("-fsanitize=bounds" CXX_COMPILER_SUPPORTS_SANITIZE_BOUNDS)
    if (CXX_COMPILER_SUPPORTS_SANITIZE_BOUNDS)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=bounds")
    endif()

    set (_saved_CRL ${CMAKE_REQUIRED_LIBRARIES})
    set (CMAKE_REQUIRED_LIBRARIES "-fsanitize=address;asan")
    CHECK_CXX_COMPILER_FLAG("-fsanitize=address" CXX_COMPILER_SUPPORTS_SANITIZE_ADDRESS)
    set (CMAKE_REQUIRED_LIBRARIES ${_saved_CRL})
    if (CXX_COMPILER_SUPPORTS_SANITIZE_ADDRESS)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    endif()

    CHECK_CXX_COMPILER_FLAG("-fsanitize=undefined" CXX_COMPILER_SUPPORTS_SANITIZE_UNDEFINED)
    if (CXX_COMPILER_SUPPORTS_SANITIZE_UNDEFINED)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined")
    endif()

    CHECK_CXX_COMPILER_FLAG("-fsanitize-address-use-after-scope" CXX_COMPILER_SUPPORTS_SANITIZE_ADDRESS_USE_AFTER_SCOPE)
    if (CXX_COMPILER_SUPPORTS_SANITIZE_ADDRESS_USE_AFTER_SCOPE)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize-address-use-after-scope")
    endif()
  endif()
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-O0")
endif() 
message("-- CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")

# column flags
set(CMAKE_CXX_FLAGS_COLUMN
    "${CMAKE_CXX_FLAGS_DEBUG}"
    CACHE STRING "Flags used by the C++ compiler during column builds."
    FORCE )

#-------------------------------------------------------------------------------
# Fortran compiler flags
#-------------------------------------------------------------------------------

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")

  set(CMAKE_Fortran_FLAGS "-g -ffree-line-length-none")

  # build type for bfb comparisons for the column physics with MPAS-Seaice
  set(CMAKE_Fortran_FLAGS_COLUMN
      "${CMAKE_Fortran_FLAGS} -O0 -m64 -fconvert=big-endian -ffree-form -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow"
      CACHE STRING "Flags used by the Fortran compiler during column builds."
      FORCE )

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")

  set(CMAKE_Fortran_FLAGS "-g -free")

else()
  message( FATAL_ERROR "-- ${CMAKE_Fortran_COMPILER_ID} fortran compiler not supported")
endif()

#-------------------------------------------------------------------------------

include(CTest)
enable_testing()
option(BUILD_TESTING "Whether to build and set up CTest tests" ON)

#RPATH settings
# CMake RPATH "always full" configuration, see:
# https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH False)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH True)
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

# NETCDF library
if (DEFINED ENV{NETCDF})
  message("-- found NETCDF=$ENV{NETCDF}")
  set(netCDF_PREFIX "$ENV{NETCDF}")
elseif (DEFINED ENV{NETCDF_ROOT})
  message("-- found NETCDF_ROOT=$ENV{NETCDF_ROOT}")
  set(netCDF_PREFIX "$ENV{NETCDF_ROOT}")
elseif (DEFINED ENV{DEMSI_NETCDF_LOC})
  message("-- found DEMSI_NETCDF_LOC=$ENV{DEMSI_NETCDF_LOC}")
  set(netCDF_PREFIX "$ENV{DEMSI_NETCDF_LOC}")
endif()
message("-- using netCDF_PREFIX=${netCDF_PREFIX}")


#find the include path and library for netCDF
find_path(netCDF_INCLUDE_DIR NAMES netcdf.h PATHS "${netCDF_PREFIX}/include")
find_library(netCDF_LIBRARY NAMES netcdf PATHS "${netCDF_PREFIX}/lib")

set(LAMMPS_PREFIX "${PROJECT_BINARY_DIR}/LAMMPS-install")
message("-- installing LAMMPS to ${LAMMPS_PREFIX}")

# determine architectures flag for LAMMPS/Kokkos 
if (DEMSI_USE_CUDA) 
  set(DEMSI_KOKKOS_ARCH "${DEMSI_GPU_ARCH};${DEMSI_CPU_ARCH}")
else()
  set(DEMSI_KOKKOS_ARCH "${DEMSI_CPU_ARCH}")
  # use the variable just so that it doesn't show up as unused at configure time
  set(DEMSI_GPU_ARCH "") 
endif()

set(LAMMPS_ARGS
 "-DCMAKE_INSTALL_PREFIX=${LAMMPS_PREFIX}"
 "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}"
 "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
 "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
 "-DPKG_MPI:BOOL=ON"
 "-DPKG_GRANULAR:BOOL=ON"
 "-DPKG_USER-DEMSI:BOOL=ON"
 "-DBUILD_LIB:BOOL=ON"
 "-DBUILD_SHARED_LIBS:BOOL=ON"
 "-DINSTALL_HEADERS:BOOL=ON"
 "-DINSTALL_LIB:BOOL=ON"
 "-DPKG_KOKKOS:BOOL=ON"
 "-DKokkos_ARCH_POWER8=ON"
 "-DKokkos_ENABLE_OPENMP:BOOL=${DEMSI_USE_OPENMP}"
 "-DKokkos_ENABLE_CUDA:STRING=${DEMSI_USE_CUDA}"
 #"-DKokkos_ARCH:STRING=${DEMSI_KOKKOS_ARCH}"
 "-DCMAKE_MACOSX_RPATH:STRING=1"
  )
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set (LAMMPS_ARGS "${LAMMPS_ARGS}" "-DCMAKE_CXX_FLAGS:STRING=-O0 -g")
  set (LAMMPS_ARGS "${LAMMPS_ARGS}" "-DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}")
  set (LAMMPS_ARGS "${LAMMPS_ARGS}" "-DKOKKOS_ENABLE_DEBUG:BOOL=ON")
else()
  set (LAMMPS_ARGS "${LAMMPS_ARGS}" "-DCMAKE_CXX_FLAGS:STRING=-O3 -g")
endif()
if (APPLE)
  set (LAMMPS_ARGS "${LAMMPS_ARGS}" "-DKOKKOS_ENABLE_LIBRT:BOOL=OFF")
endif()

message("LAMMPS_ARGS: ${LAMMPS_ARGS}")

externalproject_add(LAMMPS
    PREFIX LAMMPS
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/LAMMPS/cmake"
    CMAKE_ARGS ${LAMMPS_ARGS}
    BUILD_ALWAYS 1
    )
set(CMAKE_INSTALL_RPATH ${LAMMPS_PREFIX}/${CMAKE_INSTALL_LIBDIR} ${CMAKE_INSTALL_RPATH})
if (NOT "${CMAKE_INSTALL_LIBDIR}" STREQUAL "lib")
    #Kokkos always installs in "lib"
    set(CMAKE_INSTALL_RPATH ${LAMMPS_PREFIX}/lib ${CMAKE_INSTALL_RPATH})
endif()

set(Icepack_PREFIX "${PROJECT_BINARY_DIR}/Icepack-install")
message("-- installing Icepack to ${Icepack_PREFIX}")
add_library(Icepack
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_aerosol.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_age.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_algae.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_atmo.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_brine.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_firstyear.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_flux.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_intfc.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_itd.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_kinds.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_mechred.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_meltpond_cesm.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_meltpond_lvl.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_meltpond_topo.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_mushy_physics.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_ocean.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_orbital.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_parameters.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_shortwave.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_therm_0layer.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_therm_bl99.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_therm_itd.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_therm_mushy.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_therm_shared.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_therm_vertical.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_tracers.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_warnings.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_zbgc.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_zbgc_shared.F90
    ${PROJECT_SOURCE_DIR}/Icepack/columnphysics/icepack_zsalinity.F90
    )
set_target_properties(Icepack PROPERTIES Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/Icepack/columnphysics/" )
set_target_properties(Icepack PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${Icepack_PREFIX}")
target_include_directories(Icepack PUBLIC "${PROJECT_BINARY_DIR}/Icepack/columnphysics/")
set(CMAKE_INSTALL_RPATH ${Icepack_PREFIX} ${CMAKE_INSTALL_RPATH})

set(petsc_PREFIX "${PROJECT_BINARY_DIR}/petsc-install")
message("-- installing petsc to ${petsc_PREFIX}")

externalproject_add(petsc
    PREFIX petsc
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/petsc/"
    CONFIGURE_COMMAND cd ${PROJECT_SOURCE_DIR}/petsc/ && ./configure --prefix=${petsc_PREFIX} --download-fblaslapack
    if (APPLE)
    PETSC_ARCH=arch-darwin-c-debug
    else()
    PETSC_ARCH=arch-linux-c-debug
    endif()
    BUILD_COMMAND cd ${PROJECT_SOURCE_DIR}/petsc/ && make PETSC_DIR=${PROJECT_BINARY_DIR}/petsc all
    INSTALL_COMMAND cd ${PROJECT_SOURCE_DIR}/petsc/ && make PETSC_DIR=${PROJECT_BINARY_DIR}/petsc install
    )
set(CMAKE_INSTALL_RPATH ${petsc_PREFIX}/lib ${CMAKE_INSTALL_RPATH})

set(permon_PREFIX "${PROJECT_BINARY_DIR}/permon")
message("-- installing permon to ${permon_PREFIX}")

externalproject_add(permon
    PREFIX permon
    DEPENDS petsc
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/permon/"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND cd ${PROJECT_SOURCE_DIR}/permon/ && export PETSC_DIR=${petsc_PREFIX} && export PERMON_DIR=${PROJECT_SOURCE_DIR}/permon/ && make
    INSTALL_COMMAND ""
    )
set(CMAKE_INSTALL_RPATH ${permon_PREFIX}/lib ${CMAKE_INSTALL_RPATH})

add_subdirectory(src)

if (BUILD_TESTING)
  set(DEMSI_PYTHON "python" CACHE FILEPATH "Path to Python interpreter")
  add_subdirectory(testcases)
endif()
