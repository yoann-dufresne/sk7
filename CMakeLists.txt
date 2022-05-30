cmake_minimum_required(VERSION 3.10)

# set the project name
project(sk7)
project(sk7 VERSION 0.1)

set(CMAKE_BUILD_TYPE Debug)

# set binary output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(src/)
add_subdirectory(tests/)