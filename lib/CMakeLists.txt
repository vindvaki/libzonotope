find_package(GMP REQUIRED)
include_directories(${GMPXX_INCLUDE_DIRS})

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

add_library(zonotope_c SHARED ../src/zonotope_c.cpp)

target_link_libraries(zonotope_c ${GMP_LIBRARIES} ${GMPXX_LIBRARIES} ${CGAL_LIBRARY})

file(GLOB ZONOTOPE_C_H ../src/zonotope_c.h )

install(FILES ${ZONOTOPE_C_H} DESTINATION include)
install(TARGETS zonotope_c DESTINATION lib)
