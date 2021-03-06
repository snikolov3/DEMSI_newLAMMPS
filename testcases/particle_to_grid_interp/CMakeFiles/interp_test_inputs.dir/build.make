# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

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
CMAKE_COMMAND = /ascldap/users/projects/ppc64le/cmake/3.6.2/bin/cmake

# The command to remove a file.
RM = /ascldap/users/projects/ppc64le/cmake/3.6.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /ascldap/users/snikolo/move/DEMSI_newLAMMPS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ascldap/users/snikolo/move/DEMSI_newLAMMPS

# Utility rule file for interp_test_inputs.

# Include the progress variables for this target.
include testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/progress.make

testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs: testcases/particle_to_grid_interp/grid.nc
testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs: testcases/particle_to_grid_interp/particles_in.nc


testcases/particle_to_grid_interp/grid.nc:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/ascldap/users/snikolo/move/DEMSI_newLAMMPS/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating grid.nc, particles_in.nc"
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/particle_to_grid_interp && python /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/particle_to_grid_interp/make_testcase.py

testcases/particle_to_grid_interp/particles_in.nc: testcases/particle_to_grid_interp/grid.nc
	@$(CMAKE_COMMAND) -E touch_nocreate testcases/particle_to_grid_interp/particles_in.nc

interp_test_inputs: testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs
interp_test_inputs: testcases/particle_to_grid_interp/grid.nc
interp_test_inputs: testcases/particle_to_grid_interp/particles_in.nc
interp_test_inputs: testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/build.make

.PHONY : interp_test_inputs

# Rule to build all files generated by this target.
testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/build: interp_test_inputs

.PHONY : testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/build

testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/clean:
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/particle_to_grid_interp && $(CMAKE_COMMAND) -P CMakeFiles/interp_test_inputs.dir/cmake_clean.cmake
.PHONY : testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/clean

testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/depend:
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ascldap/users/snikolo/move/DEMSI_newLAMMPS /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/particle_to_grid_interp /ascldap/users/snikolo/move/DEMSI_newLAMMPS /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/particle_to_grid_interp /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : testcases/particle_to_grid_interp/CMakeFiles/interp_test_inputs.dir/depend

