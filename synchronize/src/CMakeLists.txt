cmake_minimum_required(VERSION 3.18)
project(Synchronize LANGUAGES CXX)

find_package(Arcane REQUIRED)

add_executable(Synchronize SynchronizeModule.cc main.cc Synchronize_axl.h)

arcane_generate_axl(Synchronize)
arcane_add_arcane_libraries_to_target(Synchronize)
target_include_directories(Synchronize PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Synchronize.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Synchronize.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Synchronize.lb.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
