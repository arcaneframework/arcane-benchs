cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

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

# ----------------------------------------------------------------------------
# Add MaHyCo and tests for MaHyCo

add_subdirectory(mahyco)
set(MAHYCO_EXE "${Mahyco_BINARY_DIR}/Mahyco")
set(MAHYCO_DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/mahyco/data")
set(MAHYCO_CPULIST 8 16 32 128 256 512 1024 2048)
foreach(X ${MAHYCO_CPULIST})
  # Used by 'mahyco/launch_test'
  set(NB_CPU ${X})
  configure_file(mahyco/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/launch_mahyco_${X}.sh @ONLY)
  add_test(NAME mahyco_${X} COMMAND /bin/sh ${CMAKE_CURRENT_BINARY_DIR}/launch_mahyco_${X}.sh)
endforeach()
