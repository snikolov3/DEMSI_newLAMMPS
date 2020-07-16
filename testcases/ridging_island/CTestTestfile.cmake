# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/ridging_island
# Build directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/ridging_island
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_ridging_island_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/smoke.py" "--testName" "smoke_ridging_island_1" "--executable" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/demsi" "--configFile" "config_ridging.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:2" "--configChanges" "remappingInterval=HOURS:1")
add_test(regression_ridging_island_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/regression.py" "--testName" "regression_ridging_island_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS" "--cmpfile" "particles_out.0001-01-01_02:00:00.nc" "--configFile" "config_ridging.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:2" "--configChanges" "remappingInterval=HOURS:1")
