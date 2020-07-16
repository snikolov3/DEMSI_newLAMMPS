# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/single_particle
# Build directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/single_particle
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_single_particle_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/smoke.py" "--testName" "smoke_single_particle_1" "--executable" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/demsi" "--configFile" "config.xml" "--nProcs" "1" "--simulationDuration" "SECONDS:10")
add_test(check_single_particle_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/check.py" "--testName" "check_single_particle_1" "--executable" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/demsi" "--configFile" "config.xml" "--nProcs" "1" "--checkscript" "check_testcase_output.py")
