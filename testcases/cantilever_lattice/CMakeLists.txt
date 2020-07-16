add_custom_command(
  OUTPUT
  grid.nc
  particles_in.nc
  COMMAND
  ${DEMSI_PYTHON}
  ${CMAKE_CURRENT_SOURCE_DIR}/make_testcase.py
  )
add_custom_target(cantilever_test_inputs ALL
  DEPENDS
  grid.nc
  particles_in.nc
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config_break.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config_break.xml
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config_nobreak.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config_nobreak.xml
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config_nobreak_equilibrium.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config_nobreak_equilibrium.xml
  COPYONLY
  )
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
endif()
add_test(
  NAME smoke_cantilever_break_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/smoke.py
  --testName smoke_cantilever_break_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config_break.xml
  --nProcs 1
  --simulationDuration HOURS:4
 )
add_test(
  NAME regression_cantilever_break_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_cantilever_break_1
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.0001-01-01_04:00:00.nc
  --configFile config_break.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration HOURS:4
  )
add_test(
  NAME smoke_cantilever_nobreak_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/smoke.py
  --testName smoke_cantilever_nobreak_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config_nobreak.xml
  --nProcs 1
  --simulationDuration HOURS:4
 )
add_test(
  NAME regression_cantilever_nobreak_1
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_cantilever_nobreak_1
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.0001-01-01_04:00:00.nc
  --configFile config_nobreak.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration HOURS:4
  )
