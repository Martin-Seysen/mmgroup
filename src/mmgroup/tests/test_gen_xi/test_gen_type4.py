"""Test C function dealing with Leech lattice vectors mod 3 of type 4

In this script we test functions dealing with vectors of type 4 in the 
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
from mmgroup.tests.chisquare import chisquare
import pytest

from mmgroup import MM0
from mmgroup.mat24 import MAT24_ORDER, ploop_theta
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.generators import gen_leech3_op_vector_atom
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.bitfunctions import binomial






#*************************************************************************
#** Convert type-4 vector mod 3 to type-4 vector mod 2
#************************************************************************/

from mmgroup.mat24 import bw24 as mat24_bw24
from mmgroup.mat24 import syndrome as mat24_syndrome
from mmgroup.mat24 import vect_to_cocode as mat24_vect_to_cocode
from mmgroup.mat24 import vect_to_gcode as mat24_vect_to_gcode
from mmgroup.mat24 import ploop_theta as MAT24_THETA_TABLE

def short_3_reduce(x):
    """Reduce coordinates of vector in Leech lattice mod 3

    The function reduces all coordinates of the vector x
    modulo 3, so that each coordinate has value 0, 1 or 2. Vector
    x must be given Leech lattice mod 3 encoding.
    """
    a = (x & (x >> 24)) & 0xffffff;
    x ^=  a | (a << 24);
    return x  & 0xffffffffffff;

def parity12(x):
    """Return parity of a 12-bit vector encoded as an integer"""
    x ^= x >> 6; x ^= x >> 3;
    return (0x96 >> (x & 7)) & 1;

def parity24(x):
    """Return parity of a 24-bit vector encoded as an integer"""
    x ^= x >> 12; x ^= x >> 6; x ^= x >> 3;
    return (0x96 >> (x & 7)) & 1;


def py_gen_leech3to2_type4(x):
    """Map type-4 vector from Leech lattice mod 3 to Leech lattice mod 2

    Here parameter ``x`` is a type-4 vector in the Leech lattice mod 3 
    in Leech lattice mod 3 encoding. 

    The function returns a type-4 vector in the Leech lattice mod 2
    corresponding to ``x`` in Leech lattice encoding.

    The result is unique. The function returns 0 if ``x`` is not of
    type 4 in the Leech lattice mod 3.

    This is a reference implementation for the C function
    **gen_leech3to2_type4** in file **gen_llech3.c**.
    """
    # uint_fast32_t gcodev, cocodev, h, w, w1, x1, syn, t, omega, res;
    x = short_3_reduce(x);
    # Let h be the support of x, i.e. the bit vector of nonzero
    # coordinates of the vector x (modulo 3)
    h = ((x >> 24) | x) & 0xffffff;
    # Let w1 and w2 be the number of indices with coordinate 1 and 2
    w = mat24_bw24(h);
    # Compute ``gcode`` and ``cocode`` for vector x. Return 0 if we 
    # detect that is not of type 4. If ``omega`` is odd then ``gcode`` 
    # has to be corrected by a term ``Omega``. At the end of the
    # switch statemnt, ``gcode`` might not correspond to a Golay
    # code vector; this means that x is not of type 4. If the scalar 
    # product of the result and ``Omega`` is one then we add a 
    # multiple of ``Omega`` to make that scalar product even.
    if w == 22:
            # type (5**1, 3**2, 1**21)
            syn = mat24_syndrome(x & 0xffffff, 0);
            gcodev = (x ^ syn) & 0xffffff;
            t = h & syn;
            cocodev = t | (0xffffff & ~h);
            if ((t == 0) or (t & (t-1))): return 0;
            omega = 0;
            #break;              
    elif w == 19:
            # type (3**5, 1**19)
            w1 = mat24_bw24(x & 0xffffff);
            x1 = x if (w1 & 1) else (x >> 24)
            x1 &= 0xffffff;
            syn = mat24_syndrome(x1, 0);
            cocodev = ~h & 0xffffff;
            if (syn & h): syn = cocodev;            
            gcodev = (x1 ^ syn) & 0xffffff;
            omega = 0;
            #break;
    elif w == 16:
            # type (2**16, 0**8)
            w1 = mat24_bw24(x & 0xffffff);
            if (w1 & 1): return 0;
            gcodev = h;
            omega = w1 >> 1;
            cocodev = x & 0xfffffff;
            #break;
    elif w in (10, 13):
            # type (4**1, 2**12, 0**11)
            # type (4**2, 2**10, 0**14)
            syn = mat24_syndrome(h & 0xffffff, 0);
            if ((h & syn) != syn): return 0;                  
            gcodev = h ^ syn;
            cocodev = syn | (x & ~syn & 0xffffff);
            w1 = mat24_bw24(cocodev);
            if (w1 & 1) : return 0;
            omega = (w1 >> 1) + parity24(syn & x) + w;
            #break; 
    elif w == 7:
            # type (6**1, 2**7, 0**16)
            syn = mat24_syndrome(h & 0xffffff, 0);
            if (syn & (syn - 1)): return 0;
            gcodev = h ^ syn;
            cocodev = (x & 0xffffff);
            w1 = mat24_bw24(cocodev);
            cocodev |=  (0 - (w1 & 1)) & syn;
            omega = ((w1 + 1) >> 1) + 1;
            #break; 
    elif w ==4:
            gcodev = 0;
            cocodev = h;
            omega = parity24(x);
            #break;
    elif w ==1:
            gcodev = cocodev = 0;
            omega = 1;
            #break;    
    else:
            return 0;        
    gcodev = mat24_vect_to_gcode(gcodev);
    if (gcodev & 0xfffff000): return 0;
    cocodev = mat24_vect_to_cocode(cocodev);
    cocodev ^= MAT24_THETA_TABLE(gcodev & 0x7ff) & 0xfff;
    # correct ``gcodev`` by term ``Omega`` if omega is odd
    gcodev ^= (omega & 1) << 11; 
    res = (gcodev << 12) ^ cocodev;
    # Correct an odd result
    if (w >= 19 and parity12(res & (res >> 12))): res ^= 0x800000;
    return res;



#####################################################################
# Creating test vectors


def rand_xsp2co1_elem(n):
    r"""Return random element of group ``G_x0`` as list of tuples

    The function returns a kind of random element of the subgroup
    ``G_x0`` of the monster as a list of tuples. That list may be
    passed to the constructor of class |MM| or |XSp2Co1|.

    Parameter ``n`` determines the number of generators 
    :math:`\xi^e, e = \pm1` in that list.
    """
    s1, l = [('p','r'),('y','r')],  [('l','n')]
    return s1 + (l + s1) * n






def create_test_elements():
    """Yield elements of the group ``G_x0`` as test data

    This generator yields elements of the group ``G_x0`` as test
    data. Each element is returned as a list o tuples as in
    function ``rand_xsp2co1_elem``. 
    """
    group_data = [
      [('x', 0x1f24), ('d', 0xf75)],
      [('x', 0x124), ('d', 0x555)],
      [('d', 0x124)],
      [('d', 0x800)],
      [('p', 187654344)],
      [('d', 0xd79), ('p', 205334671)],
      [('p', 205334671), ('d', 0xd79)],
      [('d', 0xd79), ('x', 0x1123)],
      [('y', 0x1d79)],
      [('y', 0x586)],
      [('l', 1)],
      [('l', 2)],
    ]
    for g in group_data:
        yield  g
    for n, imax  in [(1,50),(2,1000), (3,1000)]:
        for i in range(1 * imax):
            yield rand_xsp2co1_elem(n)


#####################################################################
# Test mapping of type 4 Leech vectors modulo 3 to vectors modulo 2
#####################################################################


def mul_v3(v3, g):
    """Multiply vector in Leech lattice mod 3 by element of ``G_x0``

    Here parameter ``v3`` is an element of the Leech lattice mod 3 in
    Leech lattice mod 3 encoding. Parameter ``g`` is an element of 
    the group ``G_x0`` encodes as an instance of class |MM| or |MM0|.

    The function returns the product ``v3 * g`` in  Leech lattice 
    mod 3 encoding.
    """
    result = gen_leech3_op_vector_word(v3, g._data, g.length)
    assert result & 0xffff000000000000 == 0, hex(result)
    return result
        
def mul_v2(v2, g):
    """Multiply vector in Leech lattice mod 2 by element of ``G_x0``

    Here parameter ``v2`` is an element of the Leech lattice mod 2 in
    Leech lattice encoding. Parameter ``g`` is an element of 
    the group ``G_x0`` endodes as an instance of class |MM| or |MM0|.

    The function returns the product ``v * g`` in  Leech lattice 
    encoding.
    """
    result = gen_leech2_op_word(v2, g._data, g.length)
    assert result & 0xfe000000 == 0, hex(result)
    return result

def v3_to_v2(v3):
    """Map type-4 vector from Leech lattice mod 3 to Leech lattice mod 2

    The operation of this function is equivalent to function
    ``py_gen_leech3to2_type4``. But this function is a wrapper for
    the C function ``gen_leech3to2_type4`` in file ``gen_leech3.c``. 
    """
    result = gen_leech3to2_type4(v3)
    assert result != 0, (str_v3(v3), weight_v3(v3), hex(result))
    return result


d_v3 = {0:0, 1:1, 0x1000000:2, 0x1000001:0}
def str_v3(v3):
    """Convert vector ``v3`` in Leech lattice mod 3 to a string.

    Here parameter ``v3`` is an element of the Leech lattice mod 3 in
    Leech lattice mod 3 encoding.    
    """
    l = [str(d_v3[(v3 >> i) & 0x1000001]) for i in range(24)]
    return "v3<%s>" % "".join(l)
    
w_v3 = {0:0, 1:1, 0x1000000:1, 0x1000001:0}
def weight_v3(v3):
    """Return weight of vector ``v3`` in Leech lattice mod 3.

    Here parameter ``v3`` is an element of the Leech lattice mod 3 in
    Leech lattice mod 3 encoding.    
    """
    return sum([w_v3[(v3 >> i) & 0x1000001] for i in range(24)])
    
    

@pytest.mark.gen_xi
def test_type4(verbose = 0):
    r"""Test conversion of type-4 vectors 

    Let \Omega be the type-4 vector in the Leech lattice corresponding
    to the standard frame. Let O_2 and O_3 be the images \Omega in the
    Leech lattice mod 2 and mod 3, respectively.

    For a set of elements g of the group ``G_x_0`` we convert O_3 * g to 
    a vector v2 in the Leech lattice mod 2  with function ``v3_to_v2`` 
    and we check that the result is equal to  O_2 * g. We use function
    ``create_test_elements`` for generating the elements g. 

    For the first few test data we also check this conversion against
    the reference implementation ``py_gen_leech3to2_type4``.
    """
    weights = defaultdict(int)
    for ntest, data in enumerate(create_test_elements()):
        g = MM0(data) 
        # Let \Omega be the standard type-4 vector in the Leech lattice
        # Let v3_st = \Omega in the Leech lattice mod 3
        v3_st = 1 << randint(0,23) 
        # Let v2_st = \Omega in the Leech lattice mod 3
        v2_st = 0x800000  # The standard type-4 vector \Omega      
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
        if  weights[w] <= 20:
             assert  v2 == py_gen_leech3to2_type4(v3)        
        if verbose or not ok:
            if not verbose:
                print("\nTEST %s" % (ntest+1))
                print("v3_start = " , str_v3(v3st))
                print("g =", g)
            print("v3 = v3_st*g =",str_v3(v3))
            print("weight =",w)
            print("v2 obtained= ", hex(v2))
            print("v2 expected= ", hex(v2_ref))
            if not ok:
                ERR = "Error in opation mod 3"
                raise ValueError(ERR)
    print("weights =", dict(weights))
    assert set(weights.keys()) == set([1,4,7,10,13,16,19,22])
    
    
    
#####################################################################
# Chisquare test of random type 4 Leech vectors modulo 3
#####################################################################




# The dictionary contains the number DATA_GEOMETRY[w] of type-4 
# vectors in the Leech lattice modulo 3 of weight w. 
# This table is obtained from :cite:`Iva99`, Lemma 4.4.1
DATA_GEOMETRY = {
   1: 48, 
   7: 759 * 8 * 2**7,
  22: binomial(24,3) * 3 * 2**12,
   4: binomial(24,4) * 2**4,
  10: 759 * binomial(16,2) * 2**9,
  13: 2576 * 12 * 2**12,
  19: binomial(24,5) * 2**12,
  16: 759 * 16 * 2**11,
}  
 
# Number of type-4 vectors in the Leech lattice
NUM_LEECH_TYPE4 = 398034000 
# .. must be equal to the sum of the values in DATA_GEOMETRY
assert sum(DATA_GEOMETRY.values()) ==  NUM_LEECH_TYPE4
# Inverse of the number of vectors in Leech lattice mod 3
I_NUMV3 = 3.0**(-24)


# Assemble a dictionary DICT_P for a chisquare test.
# That dictionary maps a bit weight w to itself, if that bit weight
# occurs as a weight of a type-4 vector in the Leech lattice mod 3 
# with sufficiently high probability. It maps w to 1 if w occurs as
# such a weight with low probability.
# Let P be a dictionary that maps y the probability that a random 
# vector v in the Leech lattice mod 3 has type 4 and that 
# DICT_P[weight(v)] is equal to y. Let P[0] be the probablity that
# such a random vector y is not of type 4.
DICT_P = defaultdict(int)
P = defaultdict(float)
for w, num in DATA_GEOMETRY.items():
    p = num * I_NUMV3
    assert 0 <= p < 1
    DICT_P[w] = 1 if p < 1.0/4000 else w
    P[DICT_P[w]] += p
P[0] = 1.0 - NUM_LEECH_TYPE4 * I_NUMV3
DICT_P[0] = 0

P_MIN = min([x for x in P.values() if x > 0])

RANDMOD3 = [0, 1, 0x1000000]
def rand_v3():
    """Return a random vector in the space GF(3)^{24}.

    This random vector is encoded in **Leech lattice mod 3 encoding**.
    """
    randomlist = choices(RANDMOD3, k=24)
    return sum((x << i for i, x in enumerate(randomlist)))
    
def rand_v3_dict(n):
    """Create n random vectors and group the according to their weight

    We create n random vectors in GF(3)^{24}. We return a dictionary
    ``d`` with the following entries:

    d[0] counts the vectors not of type 4 in the Leech lattice mod 3.
    d[w] counts the vectors v of type 4 such that  DICT_P[weight(v)]
    is equal to w.

    Then the value d[i] should be about P[i] / n.    
    """
    d = defaultdict(int)
    for i in range(n):
        v3 = rand_v3()
        v2 = gen_leech3to2_type4(v3) 
        if (v2 == 0):
            d[0] += 1
        else: 
            w = mat24_bw24((v3 | (v3 >> 24)) & 0xffffff)
            d[DICT_P[w]] += 1
    return d



def chisquare_v3(obtained_dict, expected_dict):
    """Perform chi square test based on dictionaries

    Dictionary ``obtained_dict`` should have the same keys as
    dictionary ``expected_dict``. The function tests the hypothesis
    that the distribution of ``obtained_dict`` is propritional to
    that of ``expected_dict`` with a chisquare test.

    The function return a pair ``(chisq, p)`` with ``chisq`` the
    chi-squared test statistic, and ``p`` the p-value of the test.
    """
    f_obt, f_exp = [],[]
    for w in obtained_dict:
        assert w in expected_dict and expected_dict[w] > 0
    factor = sum(obtained_dict.values())/sum(expected_dict.values())
    for w in expected_dict:
        if expected_dict[w] > 0:
            f_exp.append(expected_dict[w] * factor)
            f_obt.append(obtained_dict[w])
    assert min(f_exp) > 9
    chisq, p = chisquare(f_obt, f_exp = f_exp)
    return chisq, p




@pytest.mark.gen_xi
def test_chisq_type4(n = 50000, verbose = 1):
    """Test distribution of weights of type-t vectors

    The function creates a list of n random vectors in GF(3)^{24}
    and groups them with respect to the property of having type 4
    in the Leech lattice and with respect to the weight as 
    described in function ``rand_v3_dict``.

    A chisquare test fails if the p-value is less than 0.01.
    We perform at most 4 chisquare test and raise ValueError
    if all of them fail.
    """
    p_min = 0.01
    print("Check distribution of type-4 vectors mod 3") 
    for i in range(4):
        d = rand_v3_dict(n)  
        chisq, p =  chisquare_v3(d, P)
        if verbose or i or p < p_min:
            print("Chisq = %.3f, p = %.4f" % (chisq, p))
        if p >= p_min: return
    raise ValueError("Chisquare test failed") 




