"""Test C functions for computing the subtype of a Leech lattice vector

In this script we test C functions that compute the subtype of a vector
in the Leech lattice mod 2. These are functions ``gen_leech2_type`` and
``gen_leech2_type2`` in file ``gen_leech.c``.

We use the terminology defined in
the document *The C interface of the mmgroup project*, 
section *Description of the mmgroup.generators extension*.

We also count the subtypes of all 0x1000000 vectors in the Leech
lattice modulo 2 and comapare the result against the results in 
:cite:`Iva99`.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from multiprocessing import Pool
import time

import numpy as np

import pytest

from mmgroup import Xsp2_Co1, PLoop, AutPL, Cocode, MM0, MM, XLeech2
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_subtype
from mmgroup.generators import gen_leech2_coarse_subtype
from mmgroup.generators import gen_leech2_type2
from mmgroup.generators import gen_leech2_count_type2
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_type_selftest
from mmgroup.clifford12 import xsp2co1_leech2_count_type2
from mmgroup.bitfunctions import binomial



#####################################################################
# Test computation of subtype Leech lattice vector
#####################################################################

# We test the computation of the subtype of a vetor in the 
# Leech lattice with function gen_leech2_subtype().

# Next we specify certain elements of the Parker loop

# Neutral element of the Parker loop
ZERO = PLoop([])
# The element \Omega of the Parker loop
OMEGA = ~ZERO
# A (preimage of a) dodecad
DODECAD = PLoop([0,4,8, 13,14,15, 17,18,19, 21,22,23])
# A (preimage of an) octad
OCTAD = PLoop([0,1,2,3,4,5,6,7])
# A (preimage of a) complmented octad
HEXADECAD = ~OCTAD


# The following list contains a sample of a vector of each subtype 
# in the Leech lattice mod 2. Such a sample is given as a triple
#   (ploop, cocode, subytpe)
# such that  Xsp2_Co1([("x", ploop), ("d", Cocode(cocode))])
# correspond to a vector of subtype ``subtype`` in the Leech
# lattice modulo 2.
TYPE_DATA = [
    (ZERO, [], 0),
    (OMEGA, [], 0x48),
    (ZERO, [3,4], 0x20),
    (ZERO, [1], 0x21),
    (OCTAD, [], 0x22),
    (OMEGA, [1], 0x31),
    (OCTAD, [7,8], 0x34),
    (OMEGA, [1,2,3], 0x33),
    (DODECAD, [], 0x36),
    (ZERO, [0,1,2,3], 0x40),
    (OCTAD, [0,1], 0x42),
    (ZERO, [1,2,3], 0x43),
    (OCTAD, [8,9], 0x44),
    (DODECAD, [0,1], 0x46),
]


COARSE_TYPE_DATA = {
    0 : [0x00],
    1 : [0x48],
    2 : [0x20, 0x40],
    3 : [0x22, 0x42, 0x44],
    4 : [0x46],
    5 : [0x21, 0x43],
    6 : [0x34],
    7 : [0x36],
    8 : [0x31, 0x33],
}

def make_coarse_type_dict():
    d = {}
    for c, lst in COARSE_TYPE_DATA.items():
        for t in lst:
            d[t] = c
    return d

COARSE_TYPE_DICT = make_coarse_type_dict()


def rand_n_elem():
    """Return random element of the group ``N_x0``

    The group ``N_x0`` has structure ``2^{1+24}.2^{11}.M_24``.

    Multiplying a vector in the Leech lattice mod 2 with an element
    of  ``N_x0`` does not change its subtype.
    """
    return Xsp2_Co1([(x,'r') for x in "dxpy"])


def xs_vector(ploop, cocode):
    """Convert entry of table TYPE_DATA to a Leech lattice vector

    Calling ``xs_vector(ploop, cocode)``, where ``ploop`` and
    ``cocode`` are taken from an entry of the list TYPE_DATA,
    returns the vector in the Leech lattice mod 2 corresponding
    to that entry in Leech lattice encoding.
    """
    return Xsp2_Co1([("x", ploop), ("d", Cocode(cocode))]).as_xsp() 



def display_leech_vector(x):
    """Display a vector in the Leech lattice mod 2

    Here parameter ``x`` is a vector in the Leech lattice mod 2
    in Leech lattice encoding.
    """
    gcode = PLoop(x >> 12)
    bl = gcode.bit_list
    print ("GCode:\n", bl)
    if len(bl) == 16:
        l = [x for x in range(24) if not x in bl]
        pos = bl[0]
    elif len(gcode):
        pos = bl[0]
    else:
        pos = 0
    cocode = Cocode(x) + gcode.theta()
    print ("Cocode:", cocode.syndrome_list(pos))


def alternative_type2(v):
    """Alternative implementation of function ``gen_leech2_type2``

    That function returns 1 if ``v`` is of type 2 and 0 otherwise. 
    """
    a = np.array([0, v], dtype = np.uint64)
    return xsp2co1_leech2_count_type2(a, 2)

def check_leech2_subtype(x, t_expected):
    """Test computation of subtype of vector in Leech lattice mod 2

    Given a vector ``x`` in the Leech lattice mod 2 in **Leech 
    lattice encoding**, the function checkse if the subtype of ``x``
    is correctly computed as the value ``t_expected``.

    Therefore it computes the subtype with function ``gen_leech2_subtype``
    in file ``gen_leech.c`` and checks it against ``subtype``. This 
    function also checks function ``gen_leech2_type2``.   
    """
    t = gen_leech2_subtype(x)
    ok = t == t_expected
    if not ok:
        print("Error: expected Leech type: %s, obtained: %s" % (
            hex(t_expected), hex(t)))
        display_leech_vector(x)
        err = "Error in computing Leech type"
        raise ValueError(err)
    expected_type2 = t_expected if (t_expected >> 4) == 2 else 0
    found_type2 = gen_leech2_type2(x)
    ok = expected_type2 == found_type2
    if not ok:
        print("Error:  x = %s, Leech type: %s" %  
           (hex(x), hex(t_expected)), expected_type2, hex(found_type2))
        err = "Function gen_leech2_subtype2 failed"
        raise ValueError(err)
    alt_found_type2 = alternative_type2(x)
    is_type2 = (t_expected >> 4) == 2
    ok = is_type2 == alt_found_type2
    if not ok:
        print("Error:  x = %s, Leech type: %s" % (
            hex(x),  hex(t_expected)), is_type2, alt_found_type2)
        err = "Function xsp2co1_leech2_count_type2 failed"
        #raise ValueError(err)
    if not is_type2:
       assert not  found_type2 
    assert gen_leech2_type(x) == t >> 4
    ct = gen_leech2_coarse_subtype(x)
    assert ct ==  COARSE_TYPE_DICT[t], (ct, t)
  

@pytest.mark.gen_xi
def test_leech2_type(verbose = 0):
    """Test computation of subtype of vector in Leech lattice mod 2

    This function takes samples of all subtypes of vectors in
    the Leech lattice mod 2. Any such sample is transformed with
    200 random elements of the group ``N_x0``; this transformation
    does not change the subtype. The we use function 
    ``check_leech2_subtype`` for checking the computation of the
    subtype of all these transformed vectors.
    """
    for ntest, (pl, cocode, vtype) in enumerate(TYPE_DATA):
        x = xs_vector(pl, cocode)
        if verbose:
            print("\nTest %d" %(ntest+1))
            print("Expected type", hex(vtype))
            display_leech_vector(x)
        check_leech2_subtype(x, vtype)
        for i in range(200):
            g =  rand_n_elem()
            x = g.xsp_conjugate(x)
            if verbose:
                display_leech_vector(x)
            check_leech2_subtype(x, vtype)


#*************************************************************************
#** Self test for computing subypes from C file
#************************************************************************/





# The dictionary contains the number TYPE_LENGTHS[t] of vectors
# of subtype t in the Leech lattice modulo 2. This table is 
# obtained from :cite:`Iva99`, Lemmas 4.4.1 and 4.6.1.
TYPE_LENGTHS = {               # Name in :cite:`Iva99`
 0x00: 1,
 0x20: binomial(24,2) * 2,     # \Lambda_2^4
 0x21: 24 * 2**11,             # \Lambda_2^3
 0x22: 759 * 2**6,             # \Lambda_2^2
 0x31: 24 * 2**11,             # \Lambda_3^5
 0x33: binomial(24,3) * 2**11, # \Lambda_3^3
 0x34: 759 * 16 * 2**7,        # \Lambda_3^4
 0x36: 2576 * 2**10,           # \Lambda_3^2
 0x40: 2 * 1771,               # \bar{\Lambda}_4^{4a}
 0x42: 759 * 2**6,             # \bar{\Lambda}_4^{6} 
 0x43: binomial(24,3) * 2**11, # \bar{\Lambda}_4^{5}
 0x44: 15 * 759 * 2**7,        # \bar{\Lambda}_4^{4b}
 0x46: 1288 * 2**11,           # \bar{\Lambda}_4^{4c}
 0x48: 1                       # \bar{\Lambda}_4^{8}
}  
assert sum(TYPE_LENGTHS.values()) ==  0x1000000


def one_selftest_leech2(data):
    """Auxiliary function for function ``test_leech2_self``

    This is a wrapper for the C function ``gen_leech2_type_selftest``
    that counts subtypes in an interval of vectors of the Leech
    lattice modulo 2.
    """
    start, n = data
    a = np.zeros(0x50, dtype = np.uint32)
    result =  gen_leech2_type_selftest(start, n, a)
    return result, a


def gen_selftest_inputs(n):
    """Auxiliary function for function ``test_leech2_self``

    Generates inputs for function ``one_selftest_leech2`` 
    for multiprocessing.
    """
    assert 0x1000000 % n == 0
    q = 0x1000000 // n
    for  i in range(n):
        yield i*q, q


@pytest.mark.gen_xi
@pytest.mark.slow
def test_leech2_self(verbose = 0):
    """Test number of entries of all subtypes

    We count the subtypes of all 0x1000000 vector in the Leech
    lattice modulo 2, an we compile a dictionary that contains 
    the number of such vectors for each subtype. The we compare
    that dictionary against the dictionary TYPE_LENGTHS.

    We use the C function ``gen_leech2_type_selftest`` for 
    accelerating that computation, and we also use 
    multiprocessing.    
    """
    NPROCESSES = 4
    with Pool(processes = NPROCESSES) as pool:
        results = pool.map(one_selftest_leech2, 
                   gen_selftest_inputs(NPROCESSES))
    pool.join()
    result = sum(x[0] for x in results)
    a = np.zeros(0x50, dtype = np.uint32)
    for x in results:
         a += x[1]
    d = {}
    for i, n in enumerate(a):
        if n:
            d[i] = n
    if verbose:
        print("Type-4 vectors: %d\nVector types:" % result)
        for t, n in d.items():
            print(" %2x: %7d" % (t, n))
    assert sum( TYPE_LENGTHS.values() ) == 2**24
    for t, n in d.items():
        assert n == TYPE_LENGTHS[t]
    N4 = sum(n for i, n in d.items() if i & 0xf0 == 0x40)
    assert result == N4






#####################################################################
# Test counting short vectors in the Leech lattice mod 2
#####################################################################


@pytest.mark.xsp2co1
def test_leech2_count(fast = 1, verbose = 0):
    f = xsp2co1_leech2_count_type2 if fast else gen_leech2_count_type2
    a_type = np.uint64 if fast else np.uint32
    if verbose:
        print("Testing gen_leech2_type2() ... ")
    a = np.array([0]+[1 << i for i in range(24)], dtype = a_type)
    t_start = time.time()
    result = f(a, 25)
    t = time.time() - t_start
    assert result == 98280, result
    if verbose:
        print("passed, %.3f ms" % (1000*t))

#####################################################################
# Test method subtype of class Xsp2Co1
#####################################################################


def subtype_testdata():
    yield XLeech2(~PLoop())
    yield XLeech2(Cocode([0,1,2,3]))
    yield XLeech2(~PLoop(list(range(8))))
    yield XLeech2(~PLoop(list(range(8))), Cocode([8,9]))
    for i in range(50):
        yield XLeech2('r', 4)


@pytest.mark.xsp2co1
def test_subtype(verbose = 0):
    OMEGA = XLeech2(~PLoop())
    if verbose:
         print("OMEGA = ", OMEGA)
    types = set()
    for v in subtype_testdata():
         g = MM('c', v)
         v2 =  Xsp2_Co1(g)
         v2_subtype = v2.subtype
         v2ref = (OMEGA * g)
         v2ref_subtype = v2ref.subtype
         if verbose:
              print("v = ", v)
              print("g = ", g)
              print("v2 = ", v2, ", subtype =", v2_subtype)
              print("v2ref = ", v2ref, ", subtype =", v2ref_subtype)
         assert v2_subtype== v2ref_subtype, (v2_subtype, v2ref_subtype)
         types.add(v2_subtype)
    assert len(types) == 6









