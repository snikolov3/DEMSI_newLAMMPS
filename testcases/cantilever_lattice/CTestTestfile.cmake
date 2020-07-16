# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/cantilever_lattice
# Build directory: /ascldap/users/snikolo/DEMSI_cpu3/DEMSI/testcases/cantilever_lattice
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_cantilever_break_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/smoke.py" "--testName" "smoke_cantilever_break_1" "--executable" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/demsi" "--configFile" "config_break.xml" "--nProcs" "1" "--simulationDuration" "HOURS:4")
add_test(regression_cantilever_break_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/regression.py" "--testName" "regression_cantilever_break_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI" "--cmpfile" "particles_out.0001-01-01_04:00:00.nc" "--configFile" "config_break.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:4")
add_test(smoke_cantilever_nobreak_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/smoke.py" "--testName" "smoke_cantilever_nobreak_1" "--executable" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/demsi" "--configFile" "config_nobreak.xml" "--nProcs" "1" "--simulationDuration" "HOURS:4")
add_test(regression_cantilever_nobreak_1 "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/utils/testing/regression.py" "--testName" "regression_cantilever_nobreak_1" "--baseline" "" "--testdir" "/ascldap/users/snikolo/DEMSI_cpu3/DEMSI" "--cmpfile" "particles_out.0001-01-01_04:00:00.nc" "--configFile" "config_nobreak.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "HOURS:4")
