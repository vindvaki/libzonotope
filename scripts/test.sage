import zonotope
from termcolor import colored

load('zonotope.sage')

def run_test( V, identifier=""):
    n, d = V.dimensions()
    z = zonotope_from_vectors(V)

    test_results = { 'zonotope': z,
                     'generators': V }

    #
    # Test the volume
    #
    volume_libz = zonotope.zonotope_volume(V)
    volume_sage = z.volume()
    test_results['Volume'] = {
        'status' : (volume_libz == volume_sage),
        'sage': volume_sage,
        'libz': volume_libz
    }

    #
    # Test the halfspaces
    #
    halfspaces_libz = list(zonotope.zonotope_halfspaces(V))
    halfspaces_sage = [tuple(h) for h in z.Hrep_generator()]
    halfspaces_sage.sort()
    halfspaces_libz.sort()

    test_results['H-rep'] = {
        'status': halfspaces_sage == halfspaces_libz,
        'sage': halfspaces_sage,
        'libz': halfspaces_libz
    }

    #
    # Test the vertex enumeration
    #
    vertices_libz = list(zonotope.zonotope_vertices(V))
    vertices_sage = [tuple(h) for h in z.Vrep_generator()]
    vertices_sage.sort()
    vertices_libz.sort()
    test_results['V-rep'] = {
        'status': vertices_sage == vertices_libz,
        'sage': vertices_sage,
        'libz': vertices_libz,
    }

    return test_results

def handle_results(test_results, write_files = False):
    V_file_str_pattern = "fail-{0}-generators.txt"
    d, n = V.dimensions()
    status = True
    for t in ['Volume', 'H-rep', 'V-rep']:
        r = test_results[t]
        if r['status'] == False:
            status = False
            print colored('FAIL(', 'red') + t + colored(')','red'),
            if write_files == True:
                V_file_str = V_file_str_pattern.format(t)
                with open(V_file_str, 'w+') as f:
                    f.write(str(V))
    if status == True:
        print colored('passed', 'green'),
    print ""

    return status

if __name__ == '__main__':
    ntests = 100
    d, n = 2, 4
    for t in xrange(1, 1+ntests):
        print "Test {0}/{1}:".format(t, ntests),

        V = random_matrix(ZZ, n, d, x=-10, y=10)
        test_results = run_test(V, t)
        tests_passed = handle_results(test_results)

        if tests_passed == False:
            break
