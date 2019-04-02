cmake_minimum_required (VERSION 3.10)

project (CVRPTW)

set(CMAKE_CXX_STANDARD 17)

include_directories("thirdparty/SFML/include")

set(SFML_STATIC_LIBRARIES TRUE)
set(SFML_DIR "thirdparty/SFML/lib/cmake/SFML")
find_package(SFML 2.5 COMPONENTS graphics audio REQUIRED)

set(TARGET_SOURCES main.cpp)
add_executable(main ${TARGET_SOURCES})
target_link_libraries(main sfml-graphics sfml-audio)