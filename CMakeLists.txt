cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(ArcaneBenchs VERSION 1.0 LANGUAGES CXX)

# Il ne faut pas compiler la ou se trouve les sources
if(CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  message(FATAL_ERROR "You can not do in-source compilation. You have to build in a directory distinct from the source directory")
endif()

enable_testing()

# ----------------------------------------------------------------------------

# Shared libraries are mandatory
set (BUILD_SHARED_LIBS TRUE)

# ----------------------------------------------------------------------------

find_package(Arcane REQUIRED)

# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# TODO: vérifier que 'framework' existe

# ----------------------------------------------------------------------------
# Add microhydro and tests for microhydro
set(MICROHYDRO_ROOT microhydro/src)
add_subdirectory(${MICROHYDRO_ROOT} microhydro)
set(MICROHYDRO_EXE "${MicroHydro_BINARY_DIR}/MicroHydro")
set(MICROHYDRO_DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/microhydro/data")

# Get test names and number of cores from JSON file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/microhydro/test_list.json" TESTLIST_JSON)
string(JSON NB_TEST LENGTH "${TESTLIST_JSON}" "list")
math(EXPR NB_TEST "${NB_TEST} - 1")
foreach(X RANGE ${NB_TEST})
  string(JSON TEST_BASE_NAME GET "${TESTLIST_JSON}" "list" ${X} "name")
  string(JSON NB_CPU GET "${TESTLIST_JSON}" "list" ${X} "nb-core")
  set(TEST_NAME "MicroHydro.${TEST_BASE_NAME}")
  message(STATUS "MicroHydro: TEST_NAME=${TEST_NAME} NB_CPU=${NB_CPU}")
  
  configure_file(microhydro/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/launch_micro_hydro_${TEST_BASE_NAME}.sh @ONLY)
endforeach()

# ----------------------------------------------------------------------------
# Add MaHyCo and tests for MaHyCo

add_subdirectory(mahyco)
set(MAHYCO_EXE "${Mahyco_BINARY_DIR}/Mahyco")
set(MAHYCO_DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/mahyco/data")
set(MAHYCO_CPULIST 8 16 32 128 256 512 1024 2048 4096 8192)
foreach(X ${MAHYCO_CPULIST})
  # Used by 'mahyco/launch_test'
  set(NB_CPU ${X})
  configure_file(mahyco/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/launch_mahyco_${X}.sh @ONLY)
  add_test(NAME mahyco_${X} COMMAND /bin/sh ${CMAKE_CURRENT_BINARY_DIR}/launch_mahyco_${X}.sh)
endforeach()

# ----------------------------------------------------------------------------
# Add Quicksilver and tests for Quicksilver

add_subdirectory(quicksilver)
set(QUICKSILVER_EXE "${Quicksilver_BINARY_DIR}/Quicksilver")

# Get test names and number of cores from JSON file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/quicksilver/test_list.json" TESTLIST_JSON)

set(QUICKSILVER_DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/quicksilver/data/tests")
string(JSON NB_TEST LENGTH "${TESTLIST_JSON}" "tests_arc")
math(EXPR NB_TEST "${NB_TEST} - 1")
foreach(X RANGE ${NB_TEST})
  string(JSON TEST_BASE_NAME GET "${TESTLIST_JSON}" "tests_arc" ${X} "name")
  string(JSON NB_PROC GET "${TESTLIST_JSON}" "tests_arc" ${X} "nb-proc")
  string(JSON NB_THREAD GET "${TESTLIST_JSON}" "tests_arc" ${X} "nb-thread")
  set(TEST_NAME "${TEST_BASE_NAME}")
  message(STATUS "Quicksilver: TEST_NAME=${TEST_NAME} NB_PROC=${NB_PROC} NB_THREAD=${NB_THREAD}")
  configure_file(quicksilver/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/launch_quicksilver_${TEST_BASE_NAME}.sh @ONLY)
endforeach()

set(QUICKSILVER_DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/quicksilver/data/Coral2")
string(JSON NB_TEST LENGTH "${TESTLIST_JSON}" "tests_coral2")
math(EXPR NB_TEST "${NB_TEST} - 1")
foreach(X RANGE ${NB_TEST})
  string(JSON TEST_BASE_NAME GET "${TESTLIST_JSON}" "tests_coral2" ${X} "name")
  string(JSON NB_PROC GET "${TESTLIST_JSON}" "tests_coral2" ${X} "nb-proc")
  string(JSON NB_THREAD GET "${TESTLIST_JSON}" "tests_coral2" ${X} "nb-thread")
  set(TEST_NAME "${TEST_BASE_NAME}")
  message(STATUS "Quicksilver: TEST_NAME=${TEST_NAME} NB_PROC=${NB_PROC} NB_THREAD=${NB_THREAD}")
  configure_file(quicksilver/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/launch_quicksilver_${TEST_BASE_NAME}.sh @ONLY)
endforeach()

set(QUICKSILVER_DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/quicksilver/data/CTS2")
string(JSON NB_TEST LENGTH "${TESTLIST_JSON}" "tests_cts2")
math(EXPR NB_TEST "${NB_TEST} - 1")
foreach(X RANGE ${NB_TEST})
  string(JSON TEST_BASE_NAME GET "${TESTLIST_JSON}" "tests_cts2" ${X} "name")
  string(JSON NB_PROC GET "${TESTLIST_JSON}" "tests_cts2" ${X} "nb-proc")
  string(JSON NB_THREAD GET "${TESTLIST_JSON}" "tests_cts2" ${X} "nb-thread")
  set(TEST_NAME "${TEST_BASE_NAME}")
  message(STATUS "Quicksilver: TEST_NAME=${TEST_NAME} NB_PROC=${NB_PROC} NB_THREAD=${NB_THREAD}")
  configure_file(quicksilver/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/launch_quicksilver_${TEST_BASE_NAME}.sh @ONLY)
endforeach()

set(QUICKSILVER_DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/quicksilver/data/qs_original")
string(JSON NB_TEST LENGTH "${TESTLIST_JSON}" "tests_qs_original")
math(EXPR NB_TEST "${NB_TEST} - 1")
foreach(X RANGE ${NB_TEST})
  string(JSON TEST_BASE_NAME GET "${TESTLIST_JSON}" "tests_qs_original" ${X} "name")
  string(JSON NB_PROC GET "${TESTLIST_JSON}" "tests_qs_original" ${X} "nb-proc")
  string(JSON NB_THREAD GET "${TESTLIST_JSON}" "tests_qs_original" ${X} "nb-thread")
  set(TEST_NAME "${TEST_BASE_NAME}")
  message(STATUS "Quicksilver: TEST_NAME=${TEST_NAME} NB_PROC=${NB_PROC} NB_THREAD=${NB_THREAD}")
  configure_file(quicksilver/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/launch_quicksilver_${TEST_BASE_NAME}.sh @ONLY)
endforeach()
