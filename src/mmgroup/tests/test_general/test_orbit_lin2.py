"""Test C functions for the union-find algorithm.

In this module we test the classes and function in module
``mmgroup.general``
"""

from collections import defaultdict
from numbers import Integral
import pytest
import numpy as np
from copy import deepcopy
from random import sample, randint

from mmgroup import Cocode, XLeech2, PLoop, Xsp2_Co1
from mmgroup import leech2_orbits_raw
from mmgroup.bitfunctions import bitweight
from mmgroup.generators import gen_ufind_lin2_map_v_gen
from mmgroup.clifford12 import bitmatrix64_mul
from mmgroup.clifford12 import bitmatrix64_inv

from mmgroup.general import Orbit_Lin2

from mmgroup.tests.test_general.bitmatrix import chk, BitMatrix
from mmgroup.tests.test_general.bitmatrix import PermBitMatrix
from mmgroup.tests.test_general.bitmatrix import RandBitmatrix
from mmgroup.tests.test_general.bitmatrix import BitMatrix_from_Orbit_Lin2
from mmgroup.tests.test_general.bitmatrix import rand_pair_orbit
from mmgroup.tests.test_general.bitmatrix import lin2_orbits_raw
from mmgroup.tests.test_general.bitmatrix import do_test_orbit_array

def map_bitmatrix(g):
    """Mapping function for class ``GAffine``

    This function is used to map an instance of class ``GAffine``
    to an affine operation in a way suitable for class 
    ``Orbit_Lin2``
    """
    return g.map_lin2()

def lrange(*args):
    return list(range(*args))



                
                    

#####################################################################
# Test data
#####################################################################


def make_affine_group(dim = 9):
    g0 = PermBitMatrix(lrange(dim), 3)
    g1 = PermBitMatrix(lrange(1,dim) + [0])
    return Orbit_Lin2(map_bitmatrix, [g0, g1])


#####################################################################
# Main test function
#####################################################################


def do_test_affine(dim = 5, n_gen = 3, aff = True, verbose = 0):
    r"""Test the union-find algorithm on the goup H

    Here the group H is acting as a permutation group of the basis
    vectors of GF(2)^8 as described in the documentation string
    of this file.

    This function tests the functions in module ``gen_ufind_lin2.c``
    and also the python class ``Orbit_Lin2`` with the action of
    the group H on the vector space GF(2)^8.
    """
    # General test stuff
    g_list = [RandBitmatrix(dim, affine = aff)
        for i in range(n_gen)]
    G0 = Orbit_Lin2(map_bitmatrix, [g for g in g_list])
    assert G0.dim == dim

    do_test_orbit_array(G0, g_list, chk = dim <= 5, verbose = verbose)
    return


@pytest.mark.general
def test_affine(verbose = 0):
    for i, dim in enumerate([4, 5]):
        if verbose:
           print("\nTest", i+1)
        do_test_affine(dim, n_gen = 2, aff = 1, 
           verbose = verbose)

    g = make_affine_group() 
    do_test_orbit_array(g, None, chk = True, verbose = 0)



@pytest.mark.slow
@pytest.mark.general
def test_affine_slow(verbose = 0):
   do_test_affine(17, n_gen = 2, aff = 1, verbose = verbose)


"""Deprecated benchmark, do not use!!!

def test_find_bit():
    print("")
    import time
    from mmgroup.generators import gen_ufind_lin2_find_next_bit_set
    def do_test_find_bit(n, bits):
        q, r = divmod(bits, 64)
        a = np.zeros(q + 1, dtype = np.uint64)
        a[q] = 1 << r
        t_start = time.process_time()
        res = gen_ufind_lin2_find_next_bit_set(n, a, 0)
        t = (time.process_time() - t_start) / n
        assert res == n * bits
        print("Bitmap search time for %d bits: %.5f us" %
            (bits, 1000000.0*t))
    do_test_find_bit(1 << 10, 1 << 24)
    do_test_find_bit(1 << 28, 1)
"""
