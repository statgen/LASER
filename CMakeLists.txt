cmake_minimum_required(VERSION 3.8)
project(LASER VERSION 2.04)

set(LAST_MODIFIED "Feb/22/2020")

configure_file(Version.h.in Version.h)

add_compile_options(-Wall -Wpedantic)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)
set(CTEST_OUTPUT_ON_FAILURE TRUE)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

