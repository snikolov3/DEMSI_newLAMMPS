# Compute paths

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was KokkosConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#Find dependencies
INCLUDE(CMakeFindDependencyMacro)

#This needs to go above the KokkosTargets in case
#the Kokkos targets depend in some way on the TPL imports
IF(NOT TARGET Kokkos::LIBDL)
ADD_LIBRARY(Kokkos::LIBDL UNKNOWN IMPORTED)
SET_TARGET_PROPERTIES(Kokkos::LIBDL PROPERTIES
IMPORTED_LOCATION /usr/lib64/libdl.so
INTERFACE_INCLUDE_DIRECTORIES /usr/include
)
ENDIF()

GET_FILENAME_COMPONENT(Kokkos_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
INCLUDE("${Kokkos_CMAKE_DIR}/KokkosTargets.cmake")
INCLUDE("${Kokkos_CMAKE_DIR}/KokkosConfigCommon.cmake")
UNSET(Kokkos_CMAKE_DIR)
