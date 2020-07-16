add_custom_command(
  OUTPUT
  grid.nc
  forcing_fixed.nc
  particles_in.nc
  COMMAND
  ${DEMSI_PYTHON}
  ${CMAKE_CURRENT_SOURCE_DIR}/make_testcase.py
  )
add_custom_target(one_dimensional_ridging_test_inputs ALL
  DEPENDS
  grid.nc
  forcing_fixed.nc
  particles_in.nc
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config.xml
  COPYONLY
  )
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
endif()
add_test(
  NAME smoke_one_dimensional_ridging_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/smoke.py
  --testName smoke_one_dimensional_ridging_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration HOURS:2
  --configChanges "remappingInterval=HOURS:1"
  )
add_test(
  NAME regression_one_dimensional_ridging_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_one_dimensional_ridging_1
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.0001-01-01_02:00:00.nc
  --configFile config.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration HOURS:2
  --configChanges "remappingInterval=HOURS:1"
  )
