cmake_minimum_required(VERSION 3.21)
project(ArcaneSamples LANGUAGES CXX)

add_subdirectory(src)

install(DIRECTORY data reference_files DESTINATION samples/qama)
