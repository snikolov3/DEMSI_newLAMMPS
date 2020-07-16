add_custom_command(
  OUTPUT 
  grid.nc
  particles_in.nc
  COMMAND
  ${DEMSI_PYTHON}
  ${CMAKE_CURRENT_SOURCE_DIR}/make_testcase.py
  )
add_custom_target(single_particle_test_inputs ALL
  DEPENDS
  grid.nc
  particles_in.nc
  )
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
endif()
add_test(
  NAME smoke_single_particle_1
  COMMAND
  ${PROJECT_BINARY_DIR}/utils/testing/smoke.py
  --testName smoke_single_particle_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config.xml
  --nProcs 1
  --simulationDuration SECONDS:10
 )
add_test(
  NAME check_single_particle_1
  COMMAND
  ${PROJECT_BINARY_DIR}/utils/testing/check.py
  --testName check_single_particle_1
  --executable ${PROJECT_BINARY_DIR}/demsi
  --configFile config.xml
  --nProcs 1
  --checkscript check_testcase_output.py
 )
