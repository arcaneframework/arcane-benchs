cmake_minimum_required(VERSION 3.21)
project(Quicksilver LANGUAGES C CXX)

set(BUILD_SHARED_LIBS TRUE)

find_package(Arcane REQUIRED)

add_executable(Quicksilver
  main.cc
  tracking_mc/NuclearData.cc
)

add_library(QS qs/QSModule.cc)
target_include_directories(QS PUBLIC .)
target_link_libraries(QS PUBLIC arcane_core)

add_library(SamplingMC sampling_mc/SamplingMCModule.cc)
target_include_directories(SamplingMC PUBLIC .)
target_link_libraries(SamplingMC PUBLIC arcane_core)

add_library(TrackingMC tracking_mc/TrackingMCModule.cc)
target_include_directories(TrackingMC PUBLIC .)
target_link_libraries(TrackingMC PUBLIC arcane_core)

add_library(RNG rng/RNGService.cc)
target_include_directories(RNG PUBLIC .)
target_link_libraries(RNG PUBLIC arcane_core)

set(list_libs QS SamplingMC TrackingMC RNG)

target_link_libraries(Quicksilver PRIVATE ${list_libs})


set_target_properties(${list_libs} Quicksilver PROPERTIES
  INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
  INSTALL_RPATH_USE_LINK_PATH TRUE
)

include(GNUInstallDirs)

install(TARGETS ${list_libs} DESTINATION lib)
install(TARGETS Quicksilver DESTINATION bin)
install(FILES Quicksilver.config DESTINATION bin)


arcane_target_add_axl(QS qs/QS)
arcane_target_add_axl(SamplingMC sampling_mc/SamplingMC)
arcane_target_add_axl(TrackingMC tracking_mc/TrackingMC)
arcane_target_add_axl(RNG rng/RNG)


arcane_add_arcane_libraries_to_target(Quicksilver)
target_compile_options(Quicksilver PUBLIC -Wpedantic)
target_include_directories(Quicksilver PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Quicksilver.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
