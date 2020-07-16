# CMake generated Testfile for 
# Source directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/single_particle_column
# Build directory: /ascldap/users/snikolo/move/DEMSI_newLAMMPS/testcases/single_particle_column
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(regression_single_particle_column_noremap "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/regression.py" "--testName" "regression_single_particle_column_noremap" "--baseline" "" "--testdir" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS" "--cmpfile" "particles_out.1990-01-31_00.nc" "--configFile" "config.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "DAYS:30")
add_test(regression_single_particle_column_remap "/ascldap/users/snikolo/move/DEMSI_newLAMMPS/utils/testing/regression.py" "--testName" "regression_single_particle_column_remap" "--baseline" "" "--testdir" "/ascldap/users/snikolo/move/DEMSI_newLAMMPS" "--cmpfile" "particles_out.1990-01-31_00.nc" "--configFile" "config_remap.xml" "--nProcs" "1" "--nThreads" "1" "--simulationDuration" "DAYS:30")
