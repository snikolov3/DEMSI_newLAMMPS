#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "LAMMPS::lmp" for configuration "Release"
set_property(TARGET LAMMPS::lmp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(LAMMPS::lmp PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/lmp"
  )

list(APPEND _IMPORT_CHECK_TARGETS LAMMPS::lmp )
list(APPEND _IMPORT_CHECK_FILES_FOR_LAMMPS::lmp "${_IMPORT_PREFIX}/bin/lmp" )

# Import target "LAMMPS::lammps" for configuration "Release"
set_property(TARGET LAMMPS::lammps APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(LAMMPS::lammps PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/liblammps.so.0"
  IMPORTED_SONAME_RELEASE "liblammps.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS LAMMPS::lammps )
list(APPEND _IMPORT_CHECK_FILES_FOR_LAMMPS::lammps "${_IMPORT_PREFIX}/lib64/liblammps.so.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
