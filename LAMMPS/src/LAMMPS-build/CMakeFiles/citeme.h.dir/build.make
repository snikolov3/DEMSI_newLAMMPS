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

# Utility rule file for citeme.h.

# Include the progress variables for this target.
include CMakeFiles/citeme.h.dir/progress.make

CMakeFiles/citeme.h: includes/lammps/citeme.h


includes/lammps/citeme.h: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/citeme.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating includes/lammps/citeme.h"
	/home/projects/ppc64le/cmake/3.12.3/bin/cmake -E copy_if_different /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/citeme.h /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build/includes/lammps/citeme.h

citeme.h: CMakeFiles/citeme.h
citeme.h: includes/lammps/citeme.h
citeme.h: CMakeFiles/citeme.h.dir/build.make

.PHONY : citeme.h

# Rule to build all files generated by this target.
CMakeFiles/citeme.h.dir/build: citeme.h

.PHONY : CMakeFiles/citeme.h.dir/build

CMakeFiles/citeme.h.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/citeme.h.dir/cmake_clean.cmake
.PHONY : CMakeFiles/citeme.h.dir/clean

CMakeFiles/citeme.h.dir/depend:
	cd /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/cmake /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/cmake /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS/src/LAMMPS-build/CMakeFiles/citeme.h.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/citeme.h.dir/depend

