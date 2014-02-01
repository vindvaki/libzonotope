from ctypes import *

libzonotope_cdll = CDLL('./libzonotope_extern.so')

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
    return len(generators[0]), len(generators)

def _flatten_to_ctypes_array(generators, T=c_long):
    """
    Flatten generators to a C-array of type T
    """
    _generators = [x for column in generators for x in column]
    _generators = (T * (n*d))(*_generators)

#
# Exported interface
#

def zonotope_volume(generators, T=c_long):
    d, n = dimensions(generators)
    _generators = _flatten_to_ctypes_array(generators)
    return cdll_zonotope_volume_fn[T](d, n, _generators)

def zonotope_halfspaces(generators, T=c_long):
    """
    Return a generator to the list of halfspaces of the zonotope
    """
    d, n = dimensions(generators)
    _generators = _flatten_to_ctypes_array(generators)
    _halfspaces = cdll_zonotope_halfspaces_fn[T](d, n, _generators)
    return ([a[i] for i in xrange(k, k+d+1)] for k in xrange(0, n*(d+1), d+1))
