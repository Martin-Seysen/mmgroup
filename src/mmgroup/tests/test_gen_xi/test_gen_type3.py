"""Test C function dealing with Leech lattice vectors mod 3 of type 3

In this script we test functions dealing with vectors of type 3 in the 
Leech lattice modulo 3. These functions are implemented in file 
gen_leech3.c and available in the extension mmgroup.generators.

We use the terminology defined in
the document *The C interface of the mmgroup project*, 
section *Description of the mmgroup.generators extension*.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral
from collections import defaultdict
from multiprocessing import Pool, TimeoutError

import numpy as np
import pytest

from mmgroup import MM0, Xsp2_Co1
from mmgroup.mat24 import MAT24_ORDER, ploop_theta
from mmgroup.mat24 import bw24 as mat24_bw24
from mmgroup.generators import gen_leech3to2
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.generators import gen_leech3_op_vector_atom
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_type_selftest
from mmgroup.generators import gen_leech2to3_abs
from mmgroup.generators import gen_leech3_neg

from mmgroup.tests.test_gen_xi.test_gen_type4 import rand_xsp2co1_elem
from mmgroup.tests.test_gen_xi.test_gen_type4 import create_test_elements
from mmgroup.tests.test_gen_xi.test_gen_type4 import mul_v3
from mmgroup.tests.test_gen_xi.test_gen_type4 import mul_v2
from mmgroup.tests.test_gen_xi.test_gen_type4 import str_v3
from mmgroup.tests.test_gen_xi.test_gen_type4 import weight_v3
from mmgroup.tests.test_gen_xi.test_gen_type4 import chisquare_v3



#####################################################################
# Creating test vectors
#####################################################################

# Define STD_TYPE3 = -(5,1,1,...,1) to be the standard type-3 
# vector in theLeech lattice.

# Here is STD_TYPE3 (mod 3) in Leech lattice mod 3 encoding
STD_TYPE3_MOD3 = 0xfffffe000001
# Here is STD_TYPE3 (mod 2) in Leech lattice encoding
STD_TYPE3_MOD2 = 0x800800







#####################################################################
# Test mapping of type 3 Leech vectors modulo 3 to vectors modulo 2
#####################################################################


def v3_to_v2(v3):
    """Map type-3 vector from Leech lattice mod 3 to Leech lattice mod 2

    Here parameter ``v3`` is a type-3 vector in the Leech lattice mod 3 
    in Leech lattice mod 3 encoding. 

    The function returns a type-3 vector in the Leech lattice mod 2
    corresponding to ``v3`` in Leech lattice encoding.

    The result is unique. The function returns 0 if ``v3`` is not of
    type 3 in the Leech lattice mod 3
    
    This function is a wrapper for the C function ``gen_leech3to2`` 
    in file ``gen_leech3.c`` for vectors of type 3. 
    """
    result = gen_leech3to2(v3)
    assert result >> 24 == 3, (str_v3(v3), weight_v3(v3), hex(result))
    return result & 0xffffff




@pytest.mark.gen_xi
def test_type3(verbose = 0):
    r"""Test conversion of type-3 vectors 

    Let STD_TYPE3 be the type-3 vector in the Leech lattice defined
    above. Let STD_TYPE3_MOD2 and STD_TYPE3_MOD3 be the images of
    STD_TYPE3 in the Leech lattice mod 2 and mod 3, respectively.

    For a set of elements g of the group ``G_x_0`` we convert 
    STD_TYPE3_MOD3  * g to  a vector v2 in the Leech lattice mod 2  
    with function ``v3_to_v2`` and we check that the result is equal 
    to  STD_TYPE3_MOD2 * g. We use function ``create_test_elements`` 
    for generating the elements g. 
    """
    weights = defaultdict(int)
    for ntest, data in enumerate(create_test_elements()):
        g = MM0(data) 
        v3_st = STD_TYPE3_MOD3 
        v2_st = STD_TYPE3_MOD2       
        if verbose:
            print("\nTEST %s" % (ntest+1))
            print("v3_start = " , str_v3(v3_st))
            print("g =", g)
        v3 = mul_v3(v3_st, g)
        w = weight_v3(v3)  
        weights[w] += 1        
        v2_ref = mul_v2(v2_st, g) & 0xffffff
        v2 = v3_to_v2(v3) 
        ok = v2 == v2_ref 
        v3_computed = gen_leech2to3_abs(v2)
        v3neg = gen_leech3_neg(v3)
        ok_computed = ok and v3_computed in [v3, v3neg]
        if verbose or not ok or not ok_computed:
            if not verbose:
                print("\nTEST %s" % (ntest+1))
                print("v3_start = " , str_v3(v3_st))
                print("g =", g)
            print("v3 = v3_st*g =",str_v3(v3))
            print("weight =",w)
            print("v2 obtained= ", hex(v2))
            print("v2 expected= ", hex(v2_ref))
            if not ok:
                ERR = "Error in operation mod 3"
                raise ValueError(ERR)
            print("v3_computed =",str_v3(v3_computed))
            if not ok_computed:
                ERR = "Error in recomputation of v3"
                raise ValueError(ERR)


    print("weights =", dict(weights))
    assert set(weights.keys()) == set([9, 12, 21, 24])
    
    
    
#####################################################################
# Chisquare test of random type 3 Leech vectors modulo 3
#####################################################################



# The dictionary contains the number DATA_GEOMETRY[w] of type-3 
# vectors in the Leech lattice modulo 3 of weight w. 
# This table is obtained from :cite:`Iva99`, Lemma 4.4.1
DATA_GEOMETRY = {
  24: 24 * 2**12, 
   9: 759 * 16 * 2**8,
  21: 2**12 * 24 * 23 * 22 // 6,
  12: 2576 * 2**11,
}  
 
# Number of type-3 vectors in the Leech lattice
NUM_LEECH_TYPE3 = 2**24 - 2**12
# .. must be equal to the sum of the values in DATA_GEOMETRY
assert sum(DATA_GEOMETRY.values()) ==  NUM_LEECH_TYPE3
# Inverse of the number of vectors in Leech lattice mod 3
I_NUMV3 = 3.0**(-24)

# Assemble a dictionary DICT_P for a chisquare test.
# That dictionary maps a bit weight w to itself, if that bit weight
# occurs as a weight of a type-3 vector in the Leech lattice mod 3 
# with sufficiently high probability. It maps w to 1 if w occurs as
# such a weight with low probability.
# Let P be a dictionary that maps y the probability that a random 
# vector v in the Leech lattice mod 3 has type 3 and that 
# DICT_P[weight(v)] is equal to y. Let P[0] be the probablity that
# such a random vector y is not of type 3.
BLOCKSIZE = 1000000  # Minimum of type-3 vectors needed for test 
DICT_P = defaultdict(int)
MIN_P = 1.0 / 40000
P = defaultdict(float)
for w, num in DATA_GEOMETRY.items():
    p = num * I_NUMV3
    assert 0 <= p < 1
    DICT_P[w] = 1 if p < MIN_P else w
    P[DICT_P[w]] += p
P[0] = 1.0 - NUM_LEECH_TYPE3 * I_NUMV3
DICT_P[0] = 0

P_MIN = min([x for x in P.values() if x > 0])

RANDMOD3 = [0, 1, 0x1000000]
def rand_v3():
    """Return a random vector in the space GF(3)^{24}.

    This random vector is encoded in **Leech lattice mod 3 encoding**.
    """
    randomlist = choices(RANDMOD3, k=24)
    return sum((x << i for i, x in enumerate(randomlist)))
    
def rand_v3_dict(n = BLOCKSIZE):
    """Create n random vectors and group the according to their weight

    We create n random vectors in GF(3)^{24}. We return a dictionary
    ``d`` with the following entries:

    d[0] counts the vectors not of type 3 in the Leech lattice mod 3.
    d[w] counts the vectors v of type 3 such that  DICT_P[weight(v)]
    is equal to w.

    Then the value d[i] should be about P[i] / n.    
    """
    d = defaultdict(int)
    for i in range(n):
        v3 = rand_v3()
        v2 = gen_leech3to2(v3) 
        if ((v2 >> 24) != 3 or v2 == 0 ):
            d[0] += 1
        else: 
            w = mat24_bw24((v3 | (v3 >> 24)) & 0xffffff)
            d[DICT_P[w]] += 1
    return d






@pytest.mark.slow
@pytest.mark.very_slow
@pytest.mark.gen_xi
def test_chisq_type3(verbose = 0):
    """Test distribution of weights of type-t vectors

    The function creates a list of n random vectors in GF(3)^{24}
    and groups them with respect to the property of having type 3
    in the Leech lattice and with respect to the weight as 
    described in function ``rand_v3_dict``.

    A chisquare test fails if the p-value is less than 0.01.
    We perform at most 4 chisquare test and raise ValueError
    if all of them fail.
    """
    p_min = 0.01
    print("Check distribution of type-3 vectors mod 3") 
    for i in range(4):
        d = rand_v3_dict()  
        chisq, p =  chisquare_v3(d, P)
        if verbose or i or p < p_min:
            print("Chisq = %.3f, p = %.4f" % (chisq, p))
        if p >= p_min: return
    raise ValueError("Chisquare test failed") 



#####################################################################
# Test that gen_leech3to2 detect vector not of type 0,...,4
#####################################################################


# Vectors in Leech lattice mod 3 that do not map to a vector of
# type <= 4.  Source: Docmentation of the C interface, Section
# "C interface for file gen_leech_reduce_mod3.c".
VECTORS_MOD3_LARGE = [
  0x10001ff, 0x111111, 0xe, 0x1f3, 0xf1, 0xeee000
]

def rand_Co1(n_xi):
   """Return Random element of Co_2 as word of generators of G_x0

   0 <= n_ni < 4 is the number of genetors "xi" in that element
   """
   g = Xsp2_Co1('r', 'N_x0')
   for i in range(n_xi):
       g *= Xsp2_Co1([('l', randint(1,2)), ('p', 'r')])
   return g.mmdata

@pytest.mark.gen_xi
def test_v3_type_gt_4(ntests = 20):
    """Ckeck that gen_leech3to2() detects vectors of type > 4"""
    FAIL = 2**64 - 1
    for v in VECTORS_MOD3_LARGE:
        for n_xi in range(4):
            for n in range(ntests):
                g = rand_Co1(n_xi)
                v1 = gen_leech3_op_vector_word(v, g,len(g))
                assert gen_leech3to2(v1) == FAIL
