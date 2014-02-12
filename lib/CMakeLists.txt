add_library(zonotope SHARED ../src/zonotope_extern.cpp)

target_link_libraries(zonotope ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})