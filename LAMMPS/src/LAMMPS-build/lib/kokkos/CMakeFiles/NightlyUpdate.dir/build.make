# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /ascldap/users/projects/ppc64le/cmake/3.12.3/bin/cmake

# The command to remove a file.
RM = /ascldap/users/projects/ppc64le/cmake/3.12.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/cmake

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build

# Utility rule file for NightlyUpdate.

# Include the progress variables for this target.
include lib/kokkos/CMakeFiles/NightlyUpdate.dir/progress.make

lib/kokkos/CMakeFiles/NightlyUpdate:
	cd /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build/lib/kokkos && /home/projects/ppc64le/cmake/3.12.3/bin/ctest -D NightlyUpdate

NightlyUpdate: lib/kokkos/CMakeFiles/NightlyUpdate
NightlyUpdate: lib/kokkos/CMakeFiles/NightlyUpdate.dir/build.make

.PHONY : NightlyUpdate

# Rule to build all files generated by this target.
lib/kokkos/CMakeFiles/NightlyUpdate.dir/build: NightlyUpdate

.PHONY : lib/kokkos/CMakeFiles/NightlyUpdate.dir/build

lib/kokkos/CMakeFiles/NightlyUpdate.dir/clean:
	cd /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build/lib/kokkos && $(CMAKE_COMMAND) -P CMakeFiles/NightlyUpdate.dir/cmake_clean.cmake
.PHONY : lib/kokkos/CMakeFiles/NightlyUpdate.dir/clean

lib/kokkos/CMakeFiles/NightlyUpdate.dir/depend:
	cd /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/cmake /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/lib/kokkos /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build/lib/kokkos /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build/lib/kokkos/CMakeFiles/NightlyUpdate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/kokkos/CMakeFiles/NightlyUpdate.dir/depend

