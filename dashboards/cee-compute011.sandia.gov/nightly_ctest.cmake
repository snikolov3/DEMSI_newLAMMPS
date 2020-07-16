cmake_minimum_required(VERSION 3.0.1)

set(CTEST_DO_SUBMIT ON)
set(CTEST_TEST_TYPE Nightly)

set(CTEST_SITE             "cee-compute011.sandia.gov")
set(CTEST_DASHBOARD_ROOT   "/ascldap/users/daibane/demsi-nightly")
set(CTEST_CMAKE_GENERATOR  "Unix Makefiles")

set(CTEST_PROJECT_NAME "DEMSI")
set(CTEST_SOURCE_NAME  src)
set(CTEST_BUILD_NAME   "cee-compute011")
set(CTEST_BINARY_NAME  bin)

set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")
endif()
if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
endif()

configure_file("${CTEST_SCRIPT_DIRECTORY}/CTestConfig.cmake"
               "${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake" COPYONLY)

set (CTEST_COMMAND "ctest -D ${CTEST_TEST_TYPE}")

find_program(CTEST_GIT_COMMAND NAMES "git")
find_program(CXX NAMES "mpicxx")
find_program(CC NAMES "mpicc")
find_program(MPIEXEC NAMES "mpirun")

ctest_start(${CTEST_TEST_TYPE})

if(CTEST_DO_SUBMIT)
  ctest_submit(FILES "${CMAKE_CURRENT_LIST_DIR}/Project.xml"
               RETURN_VALUE SUBMIT_ERR)
  if(SUBMIT_ERR)
    message(WARNING "Cannot submit Project.xml!")
  endif()
endif()

include("${CTEST_SCRIPT_DIRECTORY}/../demsi_ctest.cmake")

demsi_ctest(
    CLEAN_BUILD
    DO_HOST
    DO_DEVICE
    BUILD_THREADS "8"
    CXX "${CXX}"
    CC "${CC}"
    NETCDF_DIR "/projects/sems/install/rhel6-x86_64/sems/tpl/netcdf/4.4.1/gcc/6.1.0/exo"
    DEMSI_PYTHON "/usr/netpub/dsbolin/Python-3.6.1/bin/python3"
    RESULT_VARIABLE ERR
    )

if (ERR)
  message(WARNING "demsi_ctest Serial returned \"${ERR}\"")
else()
  message("demsi_ctest Serial returned \"${ERR}\"")
endif()
