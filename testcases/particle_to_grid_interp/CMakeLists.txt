
# RPATH must be appended before the executable requiring that entry from RPATH is declared
set(CMAKE_INSTALL_RPATH ${PROJECT_BINARY_DIR} ${CMAKE_INSTALL_RPATH})
set(CMAKE_INSTALL_RPATH ${PROJECT_BINARY_DIR}/src ${CMAKE_INSTALL_RPATH})
add_executable(part_to_grid_test
   part_to_grid_test.cpp)

target_link_libraries(part_to_grid_test PUBLIC demsilib Icepack)

add_custom_command(
  OUTPUT
  grid.nc
  particles_in.nc
  COMMAND
  ${DEMSI_PYTHON}
  ${CMAKE_CURRENT_SOURCE_DIR}/make_testcase.py
)
add_custom_target(interp_test_inputs ALL
  DEPENDS
  grid.nc
  particles_in.nc
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.xml
  ${CMAKE_CURRENT_BINARY_DIR}/config.xml
  COPYONLY
)
add_test(
  NAME part_to_grid_interp
  COMMAND
  ${CMAKE_CURRENT_BINARY_DIR}/part_to_grid_test
)

