cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(SimpleEulerianCFD LANGUAGES C CXX)

find_package(Arcane REQUIRED)
add_executable(SimpleEulerianCFD_Original SimpleEulerianCFD_Original.cc)
arcane_add_arcane_libraries_to_target(SimpleEulerianCFD_Original)
target_compile_options(SimpleEulerianCFD_Original PRIVATE "-Wall")

# Recopie le jeu de données
configure_file(cartesian1.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
add_test(NAME simple_eulerian_cfd_original COMMAND SimpleEulerianCFD_Original cartesian1.arc)
