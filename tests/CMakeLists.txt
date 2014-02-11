find_package(GMP REQUIRED)
include_directories(${GMPXX_INCLUDE_DIRS})

add_executable(main_zonotope_halfspaces main_zonotope_halfspaces.cpp)

target_link_libraries(main_zonotope_halfspaces ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})