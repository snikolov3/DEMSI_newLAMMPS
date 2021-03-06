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

# Utility rule file for two_particles_test_inputs.

# Include the progress variables for this target.
include testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/progress.make

testcases/two_particles/CMakeFiles/two_particles_test_inputs: testcases/two_particles/grid.nc
testcases/two_particles/CMakeFiles/two_particles_test_inputs: testcases/two_particles/particles_in_bonded.nc
testcases/two_particles/CMakeFiles/two_particles_test_inputs: testcases/two_particles/particles_in_unbonded.nc


testcases/two_particles/grid.nc:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/ascldap/users/snikolo/move/DEMSI_newLAMMPS/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating grid.nc, particles_in_bonded.nc, particles_in_unbonded.nc"
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/two_particles && python /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/two_particles/make_testcase.py

testcases/two_particles/particles_in_bonded.nc: testcases/two_particles/grid.nc
	@$(CMAKE_COMMAND) -E touch_nocreate testcases/two_particles/particles_in_bonded.nc

testcases/two_particles/particles_in_unbonded.nc: testcases/two_particles/grid.nc
	@$(CMAKE_COMMAND) -E touch_nocreate testcases/two_particles/particles_in_unbonded.nc

two_particles_test_inputs: testcases/two_particles/CMakeFiles/two_particles_test_inputs
two_particles_test_inputs: testcases/two_particles/grid.nc
two_particles_test_inputs: testcases/two_particles/particles_in_bonded.nc
two_particles_test_inputs: testcases/two_particles/particles_in_unbonded.nc
two_particles_test_inputs: testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/build.make

.PHONY : two_particles_test_inputs

# Rule to build all files generated by this target.
testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/build: two_particles_test_inputs

.PHONY : testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/build

testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/clean:
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/two_particles && $(CMAKE_COMMAND) -P CMakeFiles/two_particles_test_inputs.dir/cmake_clean.cmake
.PHONY : testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/clean

testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/depend:
	cd /ascldap/users/snikolo/move/DEMSI_newLAMMPS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ascldap/users/snikolo/move/DEMSI_newLAMMPS /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/two_particles /ascldap/users/snikolo/move/DEMSI_newLAMMPS /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/two_particles /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : testcases/two_particles/CMakeFiles/two_particles_test_inputs.dir/depend

