# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/ridging_island
# Build directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/ridging_island
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_ridging_island_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/smoke.py" "--testName" "smoke_ridging_island_1" "--executable" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/demsi" "--configFile" "config_ridging.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:2" "--configChanges" "remappingInterval=HOURS:1")
add_test(regression_ridging_island_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/regression.py" "--testName" "regression_ridging_island_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI" "--cmpfile" "particles_out.0001-01-01_02:00:00.nc" "--configFile" "config_ridging.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:2" "--configChanges" "remappingInterval=HOURS:1")
