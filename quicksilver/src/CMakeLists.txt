cmake_minimum_required(VERSION 3.16)
project(Quicksilver LANGUAGES CXX)

set(BUILD_SHARED_LIBS TRUE)

find_package(Arcane REQUIRED)

add_executable(Quicksilver
  main.cc
  MC_RNG_State.cc
  tracking_mc/NuclearData.cc)


add_library(QS qs/QSModule.cc)
target_include_directories(QS PUBLIC .)
target_link_libraries(QS PUBLIC arcane_core)

add_library(SamplingMC sampling_mc/SamplingMCModule.cc)
target_include_directories(SamplingMC PUBLIC .)
target_link_libraries(SamplingMC PUBLIC arcane_core)

add_library(TrackingMC tracking_mc/TrackingMCModule.cc)
target_include_directories(TrackingMC PUBLIC .)
target_link_libraries(TrackingMC PUBLIC arcane_core)

add_library(CsvOutput simple_table_output/CsvOutputService.cc)
target_include_directories(CsvOutput PUBLIC .)
target_link_libraries(CsvOutput PUBLIC arcane_core)


target_link_libraries(Quicksilver PRIVATE QS SamplingMC TrackingMC CsvOutput)


arcane_target_add_axl(QS qs/QS)
arcane_target_add_axl(SamplingMC sampling_mc/SamplingMC)
arcane_target_add_axl(TrackingMC tracking_mc/TrackingMC)
arcane_target_add_axl(CsvOutput simple_table_output/CsvOutput)


arcane_add_arcane_libraries_to_target(Quicksilver)
target_compile_options(Quicksilver PUBLIC -Wpedantic)
target_include_directories(Quicksilver PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Quicksilver.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)