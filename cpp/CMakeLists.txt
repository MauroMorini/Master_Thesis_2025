cmake_minimum_required(VERSION 3.14)
project(cpp LANGUAGES CXX)

file(GLOB_RECURSE LIB_SOURCES src/*.cpp)
file(GLOB TEST_SOURCES tests/*.cpp)

add_library(dgfem ${LIB_SOURCES})
target_include_directories(dgfem PUBLIC 
    include/ 
    third_party/
)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CTest)
enable_testing()

add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/Catch2 EXCLUDE_FROM_ALL)
add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests PRIVATE dgfem Catch2::Catch2WithMain)

add_test(NAME fem_tests COMMAND tests)
