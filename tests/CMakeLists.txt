find_package(GMP REQUIRED)
include_directories(${GMPXX_INCLUDE_DIRS})

# test halfspaces
add_executable(test_zonotope_halfspaces test_zonotope_halfspaces.cpp)
target_link_libraries(test_zonotope_halfspaces ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

# test volume 
add_executable(test_zonotope_volume test_zonotope_volume.cpp)
target_link_libraries(test_zonotope_volume ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})