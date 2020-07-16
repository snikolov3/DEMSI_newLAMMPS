if(DEMSI_CTEST_CMAKE)
  return()
endif()
set(DEMSI_CTEST_CMAKE true)

include(${CMAKE_CURRENT_LIST_DIR}/snl_helpers.cmake)

function(demsi_ctest)
  set(BOOL_OPTS
      "CLEAN_SOURCE"
      "CLEAN_BUILD"
      "DO_HOST"
      "DO_DEVICE"
      "CHECK_FPE"
     )
  set(UNARY_OPTS
      "BUILD_THREADS"
      "RESULT_VARIABLE"
      "NETCDF_DIR"
      "CC"
      "CXX"
      "DEMSI_PYTHON"
    )
  cmake_parse_arguments(ARG "${BOOL_OPTS}" "${UNARY_OPTS}" "" ${ARGN}) 
  if (ARG_UNPARSED_ARGUMENTS)
    message(WARNING
        "demsi_ctest called with unrecognized arguments ${ARG_UNPARSED_ARGUMENTS}")
  endif()
  set(CONFIG_OPTS
      "-DCMAKE_CXX_COMPILER=${CXX}"
      "-DCMAKE_C_COMPILER=${CC}"
      )
  if (ARG_NETCDF_DIR)
    set(CONFIG_OPTS ${CONFIG_OPTS} "-DnetCDF_PREFIX:PATH=${ARG_NETCDF_DIR}")
  endif()
  if (ARG_DEMSI_PYTHON)
    set(CONFIG_OPTS ${CONFIG_OPTS} "-DDEMSI_PYTHON:FILEPATH=${ARG_DEMSI_PYTHON}")
  endif()
  set(SNL_BOOL_OPTS)
  if (ARG_CLEAN_SOURCE)
    set(SNL_BOOL_OPTS ${SNL_BOOL_OPTS} "CLEAN_SOURCE")
  endif()
  if (ARG_CLEAN_BUILD)
    set(SNL_BOOL_OPTS ${SNL_BOOL_OPTS} "CLEAN_BUILD")
  endif()
  if (ARG_DO_HOST)
    set(SNL_BOOL_OPTS ${SNL_BOOL_OPTS} "DO_UPDATE" "DO_CONFIG" "DO_BUILD")
  endif()
  if (ARG_DO_DEVICE)
    set(SNL_BOOL_OPTS ${SNL_BOOL_OPTS} "DO_TEST")
  endif()
  snl_do_subproject(${SNL_BOOL_OPTS}
      SUBPROJECT DEMSI
      REPO_URL "git@github.com:ACME-Climate/DEMSI.git"
      BRANCH "master"
      SOURCE_DIR "${CTEST_SOURCE_DIRECTORY}/DEMSI"
      BUILD_DIR "${CTEST_BINARY_DIRECTORY}/build/DEMSI"
      CONFIG_OPTS "${CONFIG_OPTS}"
      BUILD_THREADS "${ARG_BUILD_THREADS}"
      RESULT_VARIABLE SUBPROJECT_ERR
      )
  if (ARG_RESULT_VARIABLE)
    set(${ARG_RESULT_VARIABLE} ${SUBPROJECT_ERR} PARENT_SCOPE)
  endif()
endfunction(demsi_ctest)
