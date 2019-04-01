cmake_minimum_required (VERSION 3.10)

project (CVRPTW)

set(CMAKE_CXX_STANDARD 17)

set(TARGET_SOURCES main.cpp)

add_executable(main ${TARGET_SOURCES})