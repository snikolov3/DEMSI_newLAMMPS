###############################################################################
# Build documentation
###############################################################################
option(BUILD_DOC "Build LAMMPS HTML documentation" OFF)

if(BUILD_DOC)
  # Sphinx 3.x requires at least Python 3.5
  if(CMAKE_VERSION VERSION_LESS 3.12)
    find_package(PythonInterp 3.5 REQUIRED)
    set(VIRTUALENV ${PYTHON_EXECUTABLE} -m virtualenv -p ${PYTHON_EXECUTABLE})
  else()
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    if(Python3_VERSION VERSION_LESS 3.5)
      message(FATAL_ERROR "Python 3.5 and up is required to build the HTML documentation")
    endif()
    set(VIRTUALENV ${Python3_EXECUTABLE} -m virtualenv -p ${Python3_EXECUTABLE})
  endif()

  file(GLOB DOC_SOURCES ${LAMMPS_DOC_DIR}/src/[^.]*.rst)

  add_custom_command(
    OUTPUT docenv
    COMMAND ${VIRTUALENV} docenv
  )

  set(DOCENV_BINARY_DIR ${CMAKE_BINARY_DIR}/docenv/bin)

  add_custom_command(
    OUTPUT requirements.txt
    DEPENDS docenv
    COMMAND ${CMAKE_COMMAND} -E copy ${LAMMPS_DOC_DIR}/utils/requirements.txt requirements.txt
    COMMAND ${DOCENV_BINARY_DIR}/pip install -r requirements.txt --upgrade
    COMMAND ${DOCENV_BINARY_DIR}/pip install --upgrade ${LAMMPS_DOC_DIR}/utils/converters
  )

  # download mathjax distribution and unpack to folder "mathjax"
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/mathjax/es5)
    file(DOWNLOAD "https://github.com/mathjax/MathJax/archive/3.0.5.tar.gz"
      "${CMAKE_CURRENT_BINARY_DIR}/mathjax.tar.gz"
      EXPECTED_MD5 5d9d3799cce77a1a95eee6be04eb68e7)
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf mathjax.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    file(GLOB MATHJAX_VERSION_DIR ${CMAKE_CURRENT_BINARY_DIR}/MathJax-*)
    execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${MATHJAX_VERSION_DIR} ${CMAKE_CURRENT_BINARY_DIR}/mathjax)
  endif()
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/_static/mathjax)
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/mathjax/es5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/html/_static/mathjax/)

  # for increased browser compatibility
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/html/_static/polyfill.js)
    file(DOWNLOAD "https://polyfill.io/v3/polyfill.min.js?features=es6"
      "${CMAKE_CURRENT_BINARY_DIR}/html/_static/polyfill.js")
  endif()

  # note, this may run in parallel with other tasks, so we must not use multiple processes here
  add_custom_command(
    OUTPUT html
    DEPENDS ${DOC_SOURCES} docenv requirements.txt
    COMMAND ${DOCENV_BINARY_DIR}/sphinx-build -b html -c ${LAMMPS_DOC_DIR}/utils/sphinx-config -d ${CMAKE_BINARY_DIR}/doctrees ${LAMMPS_DOC_DIR}/src html
    COMMAND ${CMAKE_COMMAND} -E create_symlink Manual.html ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
  )

  # copy selected image files to html output tree
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/html/JPG)
  set(HTML_EXTRA_IMAGES balance_nonuniform.jpg balance_rcb.jpg
    balance_uniform.jpg bow_tutorial_01.png bow_tutorial_02.png
    bow_tutorial_03.png bow_tutorial_04.png bow_tutorial_05.png
    dump1.jpg dump2.jpg examples_mdpd.gif gran_funnel.png gran_mixer.png
    hop1.jpg hop2.jpg saed_ewald_intersect.jpg saed_mesh.jpg
    screenshot_atomeye.jpg screenshot_gl.jpg screenshot_pymol.jpg
    screenshot_vmd.jpg sinusoid.jpg xrd_mesh.jpg)
  set(HTML_IMAGE_TARGETS "")
  foreach(_IMG ${HTML_EXTRA_IMAGES})
    string(PREPEND _IMG JPG/)
    list(APPEND HTML_IMAGE_TARGETS "${CMAKE_CURRENT_BINARY_DIR}/html/${_IMG}")
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/${_IMG}
      DEPENDS ${LAMMPS_DOC_DIR}/src/${_IMG} ${CMAKE_CURRENT_BINARY_DIR}/html/JPG
      COMMAND ${CMAKE_COMMAND} -E copy ${LAMMPS_DOC_DIR}/src/${_IMG} ${CMAKE_BINARY_DIR}/html/${_IMG}
    )
  endforeach()

  add_custom_target(
    doc ALL
    DEPENDS html ${CMAKE_CURRENT_BINARY_DIR}/html/_static/mathjax/es5 ${HTML_IMAGE_TARGETS}
    SOURCES ${LAMMPS_DOC_DIR}/utils/requirements.txt ${DOC_SOURCES}
  )

  install(DIRECTORY ${CMAKE_BINARY_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()
