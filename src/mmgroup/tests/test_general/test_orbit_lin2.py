"""Test C functions for the union-find algorithm.

In this module we test the classe and functin in module
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
from mmgroup.clifford12 import bitmatrix64_mul
from mmgroup.clifford12 import bitmatrix64_inv

from mmgroup.general import Orbit_Lin2

from mmgroup.tests.test_general.bitmatrix import chk, BitMatrix



def PermBitMatrix(perm, b = 0):
    return BitMatrix([1 << i for i in perm], b)

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

DIM = 9

def make_affine_group():
    g0 = PermBitMatrix(lrange(DIM), 3)
    g1 = PermBitMatrix(lrange(1,DIM) + [0])
    return Orbit_Lin2(map_bitmatrix, [g0, g1])


#####################################################################
# Main test function
#####################################################################


@pytest.mark.general
def test_affine(verbose = 0):
    r"""Test the union-find algorithm on the goup H

    Here the group H is acting as a permutation group of the basis
    vectors of GF(2)^8 as described in the documentation string
    of this file.

    This function tests the functions in module ``gen_ufind_lin2.c``
    and also the python class ``Orbit_Lin2`` with the action of
    the group H on the vector space GF(2)^8.
    """
    G = make_affine_group()


