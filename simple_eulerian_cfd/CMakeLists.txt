cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(DirectCartesian LANGUAGES C CXX)

find_package(Arcane REQUIRED)
add_executable(SimpleEulerianCFD main.cc)
arcane_add_arcane_libraries_to_target(SimpleEulerianCFD)
target_compile_options(SimpleEulerianCFD PRIVATE "-Wall")

# Recopie le jeu de données
configure_file(cartesian1.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
add_test(NAME simple_eulerian_cfd COMMAND SimpleEulerianCFD cartesian1.arc)
