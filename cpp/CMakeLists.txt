cmake_minimum_required(VERSION 3.14)
project(cpp LANGUAGES CXX)

file(GLOB LIB_SOURCES src/*.cpp)

add_library(dgfem ${LIB_SOURCES})
target_include_directories(dgfem PUBLIC include/ third_party/Eigen)

include(CTest)
enable_testing()

add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/Catch2)