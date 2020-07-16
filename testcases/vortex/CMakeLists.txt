add_custom_command(
  OUTPUT
  grid.nc
  forcing_fixed.nc
  forcing_varying.0001.nc
  particles_in.nc
  particles_in_stability.nc
  COMMAND
  ${DEMSI_PYTHON}
  ${CMAKE_CURRENT_SOURCE_DIR}/make_testcase.py
  )
add_custom_target(vortex_test_inputs ALL
  DEPENDS
  grid.nc
  forcing_fixed.nc
  forcing_varying.0001.nc
  particles_in.nc
  particles_in_stability.nc
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config_fixed_forcing.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config_fixed_forcing.xml
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config_varying_forcing.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config_varying_forcing.xml
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config_fixed_forcing_stability.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config_fixed_forcing_stability.xml
  COPYONLY
  )
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
endif()
add_test(
  NAME smoke_vortex_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/smoke.py
  --testName smoke_vortex_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config_fixed_forcing.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration HOURS:1
 )
add_test(
  NAME smoke_vortex_varying_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/smoke.py
  --testName smoke_vortex_varying_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config_varying_forcing.xml
  --nProcs 1
  --nThreads 4
  --simulationDuration HOURS:1
 )
add_test(
  NAME smoke_vortex_4
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/smoke.py
  --testName smoke_vortex_4
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config_fixed_forcing.xml
  --nProcs 4
  --nThreads 2
  --simulationDuration HOURS:12
 )
add_test(
  NAME smoke_vortex_stability_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/smoke.py
  --testName smoke_vortex_stability_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config_fixed_forcing_stability.xml
  --nProcs 1
  --simulationDuration HOURS:1
 )
add_test(
  NAME regression_vortex_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_vortex_1
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.0001-01-01_04.nc
  --configFile config_fixed_forcing.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration HOURS:4
  )
add_test(
  NAME regression_vortex_varying_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_vortex_varying_1
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.0001-01-01_04.nc
  --configFile config_varying_forcing.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration HOURS:4
  )
add_test(
  NAME regression_vortex_fixed_random_4
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_vortex_fixed_random_4
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.0001-01-01_12.nc
  --configFile config_fixed_forcing_random.xml
  --nProcs 4
  --nThreads 1
  --simulationDuration HOURS:12
  )
