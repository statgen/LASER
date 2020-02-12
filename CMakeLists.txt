cmake_minimum_required(VERSION 3.8)
project(LASER VERSION 2.04)
add_compile_options(-Wall -Wpedantic)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)
set(CTEST_OUTPUT_ON_FAILURE TRUE)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

