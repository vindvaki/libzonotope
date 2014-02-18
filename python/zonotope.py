import ctypes as _ctypes

_libzonotope_cdll = _ctypes.CDLL('libzonotope_c.so')

#
# Lookup tables for dynamic typing
#

_cdll_zonotope_volume_fn = {
    _ctypes.c_long: _libzonotope_cdll.zonotope_volume_long
}

_cdll_zonotope_halfspaces_fn = {
    _ctypes.c_long: _libzonotope_cdll.zonotope_halfspaces_long
}

_cdll_zonotope_vertices_fn = {
    _ctypes.c_long: _libzonotope_cdll.zonotope_vertices_long
}

#
# Utility functions
#

def _dimensions(generators):
    if generators.dimensions is None:
        return len(generators[0]), len(generators)
    else:
        n, d = generators.dimensions()
        return d, n

def _flatten_to_ctypes_array(generators, T=_ctypes.c_long):
    """
    Flatten generators to a _ctypes array of type T
    """
    d, n = _dimensions(generators)
    generators_arr = [x for column in generators for x in column]
    generators_arr = (T * (n*d))(*generators_arr)
    return generators_arr

def _array_to_2d_list(d, n, arr):
    """
    Converts a row-major 2d-array to a 2d-list of the rows.
    """
    return [tuple(arr[i] for i in xrange(k, k+d)) for k in xrange(0, n*d, d)]

#
# Exported interface
#

def zonotope_volume(generators, T=_ctypes.c_long):
    """
    Return the volume of the zonotope
    """
    d, n = _dimensions(generators)
    generators_arr = _flatten_to_ctypes_array(generators, T)
    return _cdll_zonotope_volume_fn[T](d, n, generators_arr)

def zonotope_halfspaces(generators, T=_ctypes.c_long):
    """
    Return the list of halfspaces of the zonotope
    """
    d, n = _dimensions(generators)
    
    # set up for _ctypes function call
    generators_arr = _flatten_to_ctypes_array(generators, T)
    halfspaces_arr = _ctypes.POINTER(T)()
    _cdll_zonotope_halfspaces_fn[T].argtypes = [_ctypes.c_int,
                                                _ctypes.c_int,
                                                _ctypes.POINTER(T),
                                                _ctypes.POINTER(_ctypes.POINTER(T))]

    # obtain and format result
    num_halfspaces = _cdll_zonotope_halfspaces_fn[T](d, n, generators_arr, _ctypes.byref(halfspaces_arr))
    halfspaces = _array_to_2d_list(d+1, num_halfspaces, halfspaces_arr)

    # clean up (to avoid memory leaks)
    _libzonotope_cdll.free(halfspaces_arr)

    return halfspaces

def zonotope_vertices(generators, T=_ctypes.c_long):
    """
    Return a list of the vertices of the zonotope
    """

    d, n = _dimensions(generators)
    generators_arr = _flatten_to_ctypes_array(generators, T)
    vertices_arr = _ctypes.POINTER(T)()

    _cdll_zonotope_vertices_fn[T].argtypes = [_ctypes.c_int,
                                              _ctypes.c_int,
                                              _ctypes.POINTER(T),
                                              _ctypes.POINTER(_ctypes.POINTER(T))]

    # obtain and format the result for python
    num_vertices = _cdll_zonotope_vertices_fn[T](d, n, generators_arr, _ctypes.byref(vertices_arr))
    vertices = _array_to_2d_list(d, num_vertices, vertices_arr)

    # clean up (to avoid memory leaks)
    _libzonotope_cdll.free(vertices_arr)

    return vertices

