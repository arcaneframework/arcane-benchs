cmake_minimum_required(VERSION 3.16)
project(Quicksilver LANGUAGES CXX)

find_package(Arcane REQUIRED)

add_executable(Quicksilver 
main.cc 
QSModule.cc QS_axl.h 
SamplingMCModule.cc SamplingMC_axl.h
TrackingMCModule.cc TrackingMC_axl.h
MC_RNG_State.cc
NuclearData.cc )

arcane_generate_axl(QS)
arcane_generate_axl(SamplingMC)
arcane_generate_axl(TrackingMC)
arcane_add_arcane_libraries_to_target(Quicksilver)
target_include_directories(Quicksilver PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Quicksilver.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
