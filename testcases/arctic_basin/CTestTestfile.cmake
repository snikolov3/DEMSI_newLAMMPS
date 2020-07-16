# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/arctic_basin
# Build directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/arctic_basin
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_arctic_basin_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/smoke.py" "--testName" "smoke_arctic_basin_1" "--executable" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/demsi" "--configFile" "config.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:2" "--configChanges" "remappingInterval=HOURS:1")
add_test(regression_arctic_basin_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/regression.py" "--testName" "regression_arctic_basin_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI" "--cmpfile" "particles_arctic_basin.2000-01-01_02:00:00.nc" "--configFile" "config.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:2" "--configChanges" "remappingInterval=HOURS:1")
