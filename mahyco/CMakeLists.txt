cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
set(BUILD_SHARED_LIBS TRUE)
project(ArcaneSamples LANGUAGES C CXX)

enable_testing()

add_subdirectory(src)
