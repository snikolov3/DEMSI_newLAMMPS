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
CMAKE_SOURCE_DIR = /ascldap/users/snikolo/move/DEMSI_newLAMMPS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ascldap/users/snikolo/move/DEMSI_newLAMMPS

# Utility rule file for vortex_test_inputs.

# Include the progress variables for this target.
include testcases/vortex/CMakeFiles/vortex_test_inputs.dir/progress.make

testcases/vortex/CMakeFiles/vortex_test_inputs: testcases/vortex/grid.nc
testcases/vortex/CMakeFiles/vortex_test_inputs: testcases/vortex/forcing_fixed.nc
testcases/vortex/CMakeFiles/vortex_test_inputs: testcases/vortex/forcing_varying.0001.nc
testcases/vortex/CMakeFiles/vortex_test_inputs: testcases/vortex/particles_in.nc
testcases/vortex/CMakeFiles/vortex_test_inputs: testcases/vortex/particles_in_stability.nc


testcases/vortex/grid.nc:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/ascldap/users/snikolo/move/DEMSI_newLAMMPS/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating grid.nc, forcing_fixed.nc, forcing_varying.0001.nc, particles_in.nc, particles_in_stability.nc"
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/vortex && python /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/vortex/make_testcase.py

testcases/vortex/forcing_fixed.nc: testcases/vortex/grid.nc
	@$(CMAKE_COMMAND) -E touch_nocreate testcases/vortex/forcing_fixed.nc

testcases/vortex/forcing_varying.0001.nc: testcases/vortex/grid.nc
	@$(CMAKE_COMMAND) -E touch_nocreate testcases/vortex/forcing_varying.0001.nc

testcases/vortex/particles_in.nc: testcases/vortex/grid.nc
	@$(CMAKE_COMMAND) -E touch_nocreate testcases/vortex/particles_in.nc

testcases/vortex/particles_in_stability.nc: testcases/vortex/grid.nc
	@$(CMAKE_COMMAND) -E touch_nocreate testcases/vortex/particles_in_stability.nc

vortex_test_inputs: testcases/vortex/CMakeFiles/vortex_test_inputs
vortex_test_inputs: testcases/vortex/grid.nc
vortex_test_inputs: testcases/vortex/forcing_fixed.nc
vortex_test_inputs: testcases/vortex/forcing_varying.0001.nc
vortex_test_inputs: testcases/vortex/particles_in.nc
vortex_test_inputs: testcases/vortex/particles_in_stability.nc
vortex_test_inputs: testcases/vortex/CMakeFiles/vortex_test_inputs.dir/build.make

.PHONY : vortex_test_inputs

# Rule to build all files generated by this target.
testcases/vortex/CMakeFiles/vortex_test_inputs.dir/build: vortex_test_inputs

.PHONY : testcases/vortex/CMakeFiles/vortex_test_inputs.dir/build

testcases/vortex/CMakeFiles/vortex_test_inputs.dir/clean:
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/vortex && $(CMAKE_COMMAND) -P CMakeFiles/vortex_test_inputs.dir/cmake_clean.cmake
.PHONY : testcases/vortex/CMakeFiles/vortex_test_inputs.dir/clean

testcases/vortex/CMakeFiles/vortex_test_inputs.dir/depend:
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ascldap/users/snikolo/move/DEMSI_newLAMMPS /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/vortex /ascldap/users/snikolo/move/DEMSI_newLAMMPS /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/vortex /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/vortex/CMakeFiles/vortex_test_inputs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : testcases/vortex/CMakeFiles/vortex_test_inputs.dir/depend

