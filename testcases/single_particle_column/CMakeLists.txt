add_custom_command(
  OUTPUT 
  forcing_CORE_monthly_clim_atmos.0000.nc
  forcing_CORE_monthly_clim_ocean.0000.nc
  forcing_CORE_six_hourly.nc
  grid.nc
  particles_in.nc
  COMMAND
  ${DEMSI_PYTHON}
  ${CMAKE_CURRENT_SOURCE_DIR}/make_testcase.py
  )
add_custom_target(single_particle_column_test_inputs ALL
  DEPENDS
  forcing_CORE_monthly_clim_atmos.0000.nc
  forcing_CORE_monthly_clim_ocean.0000.nc
  forcing_CORE_six_hourly.nc
  grid.nc
  particles_in.nc
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config.xml
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config_remap.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config_remap.xml
  COPYONLY
  )
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output_noremap)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_noremap)
endif()
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output_remap)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_remap)
endif()
add_test(
  NAME regression_single_particle_column_noremap
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_single_particle_column_noremap
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.1990-01-31_00.nc
  --configFile config.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration DAYS:30
  )
add_test(
  NAME regression_single_particle_column_remap
  COMMAND
  ${PROJECT_SOURCE_DIR}/utils/testing/regression.py
  --testName regression_single_particle_column_remap
  --baseline "$ENV{DEMSI_TEST_BASELINE}"
  --testdir ${PROJECT_BINARY_DIR}
  --cmpfile particles_out.1990-01-31_00.nc
  --configFile config_remap.xml
  --nProcs 1
  --nThreads 1
  --simulationDuration DAYS:30
  )
