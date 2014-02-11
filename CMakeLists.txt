cmake_minimum_required(VERSION 2.8)

project (LIBZONOTOPE)

set(CMAKE_CXX_FLAGS "-std=c++11")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_subdirectory(tests)
add_subdirectory(src)