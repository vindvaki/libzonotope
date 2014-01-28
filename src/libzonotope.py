from ctypes import *

libzonotope_cdll = CDLL('./libzonotope.so')

def zonotope_volume(generators):
    d = len(generators[0])
    n = len(generators)
    _generators = [x for column in generators for x in column]
    _generators = (c_long * (n*d))(*_generators)
    return libzonotope_cdll._zonotope_volume(d,n, _generators)
