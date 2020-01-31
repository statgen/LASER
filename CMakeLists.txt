cmake_minimum_required(VERSION 3.8)
project(LASER VERSION 2.04)
add_compile_options(-Wall -Wpedantic)

set(CMAKE_VERBOSE_MAKEFILE ON)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

