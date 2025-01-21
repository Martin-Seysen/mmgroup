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
from mmgroup.tests.general.test_gen_ufind import chk

class GAffine:
    """Model an affine operation on a vector space.

    The constructor takes two parameters ``a``, ``b``, with ``a``
    :math:`n \times n` bit matrix given as an array of ``n`` bit
    vectors, and ``b`` a bit vector. Bit vectors are encoded as
    integers. Such a pair encodes the mapping
    :math:`v \mapsto v \cdot a + b`.

    Instances of this class may be multiplied or raised to the
    power 0, 1, or -1. They may be left multipled with an
    integer denoting a bit vector.

    Member ``dim`` contains the dimension of the matrix ``a``.
    """

    def __init__(self, a, b = 0):
        self.dim = len(a)
        assert 1 <= self.dim <= 24
        self.m = np.array(a, dtype = np.uint32)
        self.m = np.append(self.m, np.uint32(b & 0xffffffff))
        self.m &= self.mask
        self.m[self.dim] |= 1 << self.dim
        
    @property 
    def mask(self):
        return (1 << self.dim) - 1

    def map_lin2(self):
        return list(self.m[:self.dim]), self.m[self.dim] & self.mask

    def __mul__(self, other):
        if not isinstance(other, GAffine):
            return NotImplemented
        assert self.dim == other.dim
        a = np.array(self.m, dtype = np.uint64)
        b = np.array(other.m, dtype = np.uint64)
        n = self.dim + 1
        bitmatrix64_mul(a, b, n, n, a)
        return GAffine(a[:self.dim], a[self.dim] & self.mask)

    def invert(self, other):
        a = np.array(self.m, dtype = np.uint64)
        bitmatrix64_inv(a, self.dim + 1)
        return GAffine(a[:self.dim], a[self.dim] & self.mask)

    def neutral(self):
        return GAffine([1 << i for i in range(self.dim)])

    def __rmul__(self, other):
        if not isinstance(other, Integral):
            return NotImplemented
        a = np.array(self.m, dtype = np.uint64)
        other |= 1 << self.dim
        res = bitmatrix64_vmul(other, a, self.dim + 1)
        return other & self.mask

    def __pwr__(self, exp):
        if exp == 1:
            return self
        elif exp == 0:
            return self.neutral()
        elif exp == -1:
            return self.invert()
        else:
            return NotImplemented


def PermAffine(perm, b = 0):
    return GAffine([1 << i for i in perm], b)

def map_gaffine(g):
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
    g0 = PermAffine(lrange(DIM), 3)
    g1 = PermAffine(lrange(1,DIM) + [0])
    return Orbit_Lin2(map_gaffine, [g0, g1])


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


