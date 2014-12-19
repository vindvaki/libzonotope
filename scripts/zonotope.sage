#!/usr/bin/env sage

"""
A collection of functions for working with zonotopes in Sage.
"""

def zonotope_from_segments(segments):
    """
    Returns the polyhedron for the zonotope given by the Minkowski sum
    of ``segments``, where ``segments`` is a list of generating
    line-segment endpoint points pairs ``(s,t)``.
    """
    return sum(Polyhedron(vertices=s) for s in segments)

def zonotope_from_vectors(V):
    """
    Returns the zonotope constructed from the vectors ``v`` iterated
    over ``V``, using the line segments ``[0,v]`` as generators.
    """
    return zonotope_from_segments(segments_from_origin(V))

def zonotope_volume(V):
    """
    Returns the volume of the zonotope with generators ``V``.
    """
    n, d = V.dimensions()
    return sum(abs(Matrix(U).det()) for U in Combinations(V, d))

def segments_from_origin(V):
    """
    Convert an ``n-by-d`` matrix ``V`` to an iterator of segments of
    the form ``(zero_vector(d),v)`` for the rows ``v`` in ``V``.
    """
    n, d = V.dimensions()
    origin = zero_vector(d)
    return ((origin, row) for row in V)

def robust_region(V, z = None):
    """
    Computes the robust sub-polytope of the zonotope ``z`` that was
    generated from the matrix ``V``.
    """
    
    if z is None:
        z = zonotope_from_vectors(V)
    
    ieqs = []
    for H in z.Hrepresentation():
        # We have H = { x : a*x + b >= 0 } so
        # 
        #   H - v = { x - v: a*x + b >= 0 }
        #         = { y : a*y + (b + a*v) >= 0 }
        #
        #  so we only need the endpoint v minimizing a*v
        a = H.A()
        b = H.b() + min(a*v for v in V)
        ieqs.append([b] + list(a))
    return Polyhedron(ieqs=ieqs)

def robust_region_naive(V, z = None):
    """
    See ``robust_region(V, z)``.
    """
    
    if z is None:
        z = zonotope_from_vectors(V)
    
    n, d = V.dimensions()
    origin = zero_vector(d)
    return reduce(lambda a,b: a&b, (z-v for v in V))


def zonotope_halfspaces(V, verbose=False):
    """
    A proof of concept implementation that generates the
    H-representation of the zonotope directly from the generators.
    """

    n, d = V.dimensions()
    ieqs = set([])
    generator_sum = sum(V)

    for U in Combinations(V, d-1):
        normal = Matrix(ZZ, U).transpose().integer_kernel().gen()
        offset_vector = zero_vector(ZZ, d)
        for v in V:
            normal_proj = normal.dot_product(v)
            if normal_proj <= 0:
                offset_vector += v
        
        # standardize the normal to the smallest integral vector representing
        # the same hyperplane
        standardizer = abs(gcd(normal))
        if standardizer > 0:
            normal /= standardizer
        
        offset = - ( normal.dot_product(offset_vector) )
        # { x : normal*x + offset >= 0 } is one of the halfspaces
        
        # comptue the antipodal inequality
        antipode_normal = -normal
        antipode_offset = ( generator_sum.dot_product(normal) ) + offset
        # { x: antipode_normal + antipode_offset >= 0 } is the antipodal halfspace
        
        ieqs.add(tuple(offset) + tuple(normal))
        ieqs.add(tuple(antipode_offset) + tuple(antipode_normal))
        
    if verbose == True:
        print "n={n}, d={d}, n_ieqs={n_ieqs}".format(\
            n=n, d=d, n_ieqs=len(ieqs))
          
    return ieqs


def data_iterator_from_file(f, number_field=QQ):
    """
    Converts a data file formatted according to make_data to an iterator over
    its numerical values.
    """
    return (number_field(i) for i in f.read().split())

def generators_from_iterator(data, number_field=QQ):
    """
    Reads a data iterator and returns an iterator over the matrices defined by
    the iterator.
    """
    n_tests = data.next()
    t = 0
    while t < n_tests:
        n = data.next()
        d = data.next()
        generators = matrix(number_field, n, d)
        for i in xrange(n):
            for j in xrange(d):
                generators[i,j] = data.next()
        yield generators
        t += 1

def zonotope_halfspaces_test(f, number_field=QQ):
    data_iterator = data_iterator_from_file(f, number_field)
    for V in generators_from_iterator(data_iterator):
        zonotope_halfspaces(V, verbose=True)

def show_robust_region(z, V, alpha_baseline=0.6, **kwargs):
    """
    Returns an object for the visualization of the generalized robust
    regions of ``z`` with respect to the generating matrix ``V``.
    """

    n, d = V.dimensions()
    assert d <= 3

    # Init opacity options
    alpha = alpha_baseline / n
    if d <= 2:
        kw_alpha = 'alpha'
    if d == 3:
        kw_alpha = 'opacity'
    kwargs[kw_alpha] = alpha

    # Handle scaling of outlines etc for resolution-independence in
    # the output.
    default_figsize = 5
    if 'figsize' not in kwargs:
        kwargs['figsize'] = default_figsize
#    fig_scale = float(kwargs['figsize']) / default_figsize
#    kwargs['thickness'] = 3 * fig_scale
    
    # Overlay all sub-zonotopes when a generator from V is removed to
    # the visualization. We limit the opacity of each sub-zonotope,
    # such intersections and outlines are easier to discern.
    output = z.show(**kwargs)
    robust = z
    origin = zero_vector(d)
    for v in V:
        z_v = z - v
        output += (z & z_v).show(**kwargs)
        robust &= z_v
    
    # The robust region should stand out, so we give it full opacity
    kwargs[kw_alpha] = 1
    output += robust.show(**kwargs)

    # Similarly, the generators need to stand out, so they get greater
    # thickness (there is only limited built in support for changing
    # the color, so we don't do that for now)
    generator_thickness = 6
#    kwargs['thickness'] = generator_thickness * fig_scale
    if d == 2:
        kwargs['color'] = 'red'
    for v in V:
        segment = Polyhedron(vertices=(origin,v))
        output += segment.show(**kwargs)

    return output

#if __name__ == '__main__':
#    # Set up a 2D demo zonotope. Works pretty much as expected.
#    V_2 = matrix(QQ, [[-1,1],[1,0], [1,-2], [2,3], [-0.3,-1]])
#    z_2 = zonotope_from_vectors(V_2)
#    r_2 = show_robust_region(z_2, V_2)
#    r_2.save('robust_2.png', axes=False, figsize=40)
#    
#    # Set up a 3D demo zonotope. This is still broken -- the robust
#    # region graphics object doesn't seem to play well with the
#    # raytracer that Sage uses.
#    V_3 = matrix(QQ, [[-1,1,1],[1,0,1], [1,-5,1], [2,3,-1], [1,-1,-1]])
#    z_3 = zonotope_from_vectors(V_3)
#    r_3 = show_robust_region(z_3, V_3, frame=False)
#    r_3.save('robust_3.png')
