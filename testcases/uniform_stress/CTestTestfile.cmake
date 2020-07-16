# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/uniform_stress
# Build directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/uniform_stress
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_uniform_stress_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/smoke.py" "--testName" "smoke_uniform_stress_1" "--executable" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/demsi" "--configFile" "config.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "SECONDS:1500")
add_test(regression_uniform_stress_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/regression.py" "--testName" "regression_uniform_stress_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI" "--cmpfile" "particles_out.0001-01-01_00:25:00.nc" "--configFile" "config.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "SECONDS:1500")
