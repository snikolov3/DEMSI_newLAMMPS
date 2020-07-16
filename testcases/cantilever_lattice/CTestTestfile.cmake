# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/cantilever_lattice
# Build directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/cantilever_lattice
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_cantilever_break_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/smoke.py" "--testName" "smoke_cantilever_break_1" "--executable" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/demsi" "--configFile" "config_break.xml" "--nProcs" "1" "--simulationDuration" "HOURS:4")
add_test(regression_cantilever_break_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/regression.py" "--testName" "regression_cantilever_break_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS" "--cmpfile" "particles_out.0001-01-01_04:00:00.nc" "--configFile" "config_break.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:4")
add_test(smoke_cantilever_nobreak_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/smoke.py" "--testName" "smoke_cantilever_nobreak_1" "--executable" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/demsi" "--configFile" "config_nobreak.xml" "--nProcs" "1" "--simulationDuration" "HOURS:4")
add_test(regression_cantilever_nobreak_1 "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/regression.py" "--testName" "regression_cantilever_nobreak_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS" "--cmpfile" "particles_out.0001-01-01_04:00:00.nc" "--configFile" "config_nobreak.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:4")
