/**
 * @brief zonotope_volume_T Compute the volume of a zonotope
 * @param d dimension of the output space
 * @param n number of generators
 * @param generators a d-by-n column major matrix of the generators
 * @return the volume of the zonotope
 */
long   zonotope_volume_long(  int d, int n, const long*   generators);
double zonotope_volume_double(int d, int n, const double* generators);

/**
 * @brief zonotope_halfspaces_T Compute the hyperplane representation of a zonotope
 * @param d the dimension of the output space
 * @param n the number of generators
 * @param generators a d-by-n column major matrix of the generators
 * @param halfspaces on exit, *halfspaces points to the list of bounding halfspaces
 * @return the number of halfspaces in the output
 */
long zonotope_halfspaces_long(  const int d, const int n, const long*   generators, long**   halfspaces);
long zonotope_halfspaces_double(const int d, const int n, const double* generators, double** halfspaces);

/**
 * @brief zonotope_vertices_T Compute the vertices of a zonotope
 * @param d the dimension of the output space
 * @param n the number of generators
 * @param generators a d-by-n column major matrix of the generators
 * @param vertices on exit, *vertices points to the list of vertices of the zonotope
 * @return the number of vertices in the output
 */
long zonotope_vertices_long(  const int d, const int n, const long*   generators, long**   vertices);
long zonotope_vertices_double(const int d, const int n, const double* generators, double** vertices);
