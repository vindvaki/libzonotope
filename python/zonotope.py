from ctypes import *
from ctypes.util import find_library

libzonotope_cdll = CDLL('libzonotope_c.so')

#
# Lookup tables for dynamic typing
#

cdll_zonotope_volume_fn = {
    c_long: libzonotope_cdll.zonotope_volume_long
}

cdll_zonotope_halfspaces_fn = {
    c_long: libzonotope_cdll.zonotope_halfspaces_long
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

def _flatten_to_ctypes_array(generators, T=c_long):
    """
    Flatten generators to a ctypes array of type T
    """
    d, n = _dimensions(generators)
    _generators = [x for column in generators for x in column]
    _generators = [x for x in _generators]
    _generators = (T * (n*d))(*_generators)
    return _generators

def _array_to_2d_list(d, n, arr):
    """
    Converts a row-major 2d-array to a 2d-list of the rows.
    """
    return [tuple(arr[i] for i in xrange(k, k+d)) for k in xrange(0, n*d, d)]

#
# Exported interface
#

def zonotope_volume(generators, T=c_long):
    d, n = _dimensions(generators)
    _generators = _flatten_to_ctypes_array(generators)
    return cdll_zonotope_volume_fn[T](d, n, _generators)

def zonotope_halfspaces(generators, T=c_long):
    """
    Return a generator to the list of halfspaces of the zonotope
    """
    d, n = _dimensions(generators)
    print d,n
    
    # set up for ctypes function call
    _generators = _flatten_to_ctypes_array(generators)
    _halfspaces = POINTER(T)()
    cdll_zonotope_halfspaces_fn[T].argtypes = [c_int, c_int, POINTER(c_long), POINTER(POINTER(T))]

    # obtain and format result
    num_halfspaces = cdll_zonotope_halfspaces_fn[T](d, n, _generators, byref(_halfspaces))
    halfspaces = _array_to_2d_list(d+1, num_halfspaces, _halfspaces)

    # clean up (to avoid memory leaks)
    libzonotope_cdll.free(_halfspaces)

    return halfspaces
