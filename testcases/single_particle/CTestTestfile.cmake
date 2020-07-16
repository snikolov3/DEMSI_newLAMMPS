# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/single_particle
# Build directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/single_particle
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_single_particle_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/smoke.py" "--testName" "smoke_single_particle_1" "--executable" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/demsi" "--configFile" "config.xml" "--nProcs" "1" "--simulationDuration" "SECONDS:10")
add_test(check_single_particle_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/check.py" "--testName" "check_single_particle_1" "--executable" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/demsi" "--configFile" "config.xml" "--nProcs" "1" "--checkscript" "check_testcase_output.py")
