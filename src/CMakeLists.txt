add_library(zonotope_extern SHARED zonotope_volume_extern.cpp)

target_link_libraries(zonotope_extern ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})