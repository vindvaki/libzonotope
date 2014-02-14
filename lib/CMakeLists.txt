add_library(zonotope_c SHARED ../src/zonotope_extern.cpp)

target_link_libraries(zonotope_c ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

install(FILES ../src/zonotope_c.h DESTINATION include)
install(TARGETS zonotope_c DESTINATION lib)