cmake_minimum_required(VERSION 3.18)
project(ArcaneLoopsBench LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(ArcaneLoopsBench ArcaneLoopsModule.cc main.cc ArcaneLoops_axl.h)

arcane_generate_axl(ArcaneLoops)
arcane_add_arcane_libraries_to_target(ArcaneLoopsBench)
target_include_directories(ArcaneLoopsBench PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(ArcaneLoopsBench.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
