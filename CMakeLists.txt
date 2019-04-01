cmake_minimum_required (VERSION 3.10)

project (CVRPTW)

set(CMAKE_CXX_STANDARD 17)

#include_directories("thirdparty/SFML/include")

set(TARGET_SOURCES main.cpp)
add_executable(main ${TARGET_SOURCES})
#target_link_libraries(main  "thirdparty/SFML/lib/sfml-graphics" "thirdparty/SFML/lib/sfml-window"  "thirdparty/SFML/lib/sfml-system")