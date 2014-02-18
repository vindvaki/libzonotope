cmake_minimum_required(VERSION 2.8)

project (LIBZONOTOPE)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include_directories(include)

file( GLOB_RECURSE LIBZONOTOPE_HPP include/*.hpp )
file( GLOB LIBZONOTOPE_PY python/*.py )
add_custom_target( libzonotope_headers SOURCES ${LIBZONOTOPE_HPP} )
add_custom_target( libzonotope_python  SOURCES ${LIBZONOTOPE_PY} )

#
# target subdirectories
#

add_subdirectory(tests)
add_subdirectory(lib)
#add_subdirectory(python)

#
# ``uninstall`` target
# 

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
               IMMEDIATE @ONLY)

add_custom_target(uninstall
                  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
