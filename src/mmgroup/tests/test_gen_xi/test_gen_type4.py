from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral
from collections import defaultdict
from multiprocessing import Pool, TimeoutError

import numpy as np
import scipy
import scipy.special
from scipy.stats import chisquare
import pytest

from mmgroup import MM
from mmgroup.mat24 import MAT24_ORDER, ploop_theta
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.generators import gen_leech3_op_vector_atom
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_type_selftest


OMEGA0_V3  = 0x1




#*************************************************************************
#** Convert type-4 vector mod 3 to type-4 vector mod 2
#************************************************************************/

from mmgroup.mat24 import bw24 as mat24_bw24
from mmgroup.mat24 import syndrome as mat24_syndrome
from mmgroup.mat24 import vect_to_cocode as mat24_vect_to_cocode
from mmgroup.mat24 import vect_to_gcode as mat24_vect_to_gcode
from mmgroup.mat24 import ploop_theta as MAT24_THETA_TABLE

def short_3_reduce(x):
    a = (x & (x >> 24)) & 0xffffff;
    x ^=  a | (a << 24);
    return x  & 0xffffffffffff;

def parity12(x):
    x ^= x >> 6; x ^= x >> 3;
    return (0x96 >> (x & 7)) & 1;

def parity24(x):
    x ^= x >> 12; x ^= x >> 6; x ^= x >> 3;
    return (0x96 >> (x & 7)) & 1;

def cond(c, a, b):
    """Substitute for the C conditional   'c ? a : b' """
    return a if c else b

def py_gen_leech3to2_type4(x):
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
            x1 = cond((w1 & 1) , x , (x >> 24)) & 0xffffff;
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
#####################################################################


def rand_xsp2co1_elem(n):
    s1, l = [('p',),('y',)],  [('l','n')]
    return s1 + (l + s1) * n






def create_test_elements():
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
    assert g in MM
    result = gen_leech3_op_vector_word(v3, g._data, g.length)
    assert result & 0xffff000000000000 == 0, hex(result)
    return result
        
def mul_v2(v2, g):
    assert g in MM
    result = gen_leech2_op_word(v2, g._data, g.length)
    assert result & 0xfe000000 == 0, hex(result)
    return result

def v3_to_v2(v3):
    result = gen_leech3to2_type4(v3)
    assert result != 0, (str_v3(v3), weight_v3(v3), hex(result))
    return result


d_v3 = {0:0, 1:1, 0x1000000:2, 0x1000001:0}
def str_v3(v3):
    l = [str(d_v3[(v3 >> i) & 0x1000001]) for i in range(24)]
    return "v3<%s>" % "".join(l)
    
w_v3 = {0:0, 1:1, 0x1000000:1, 0x1000001:0}
def weight_v3(v3):
    return sum([w_v3[(v3 >> i) & 0x1000001] for i in range(24)])
    
    

@pytest.mark.gen_xi
def test_type4(verbose = 0):
    weights = defaultdict(int)
    for ntest, data in enumerate(create_test_elements()):
        g = MM(*data) 
        v3_st = 1 << randint(0,23) 
        v2_st = 0x800000        
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



def binom(n, k):
    return int(scipy.special.binom(n, k) + 0.1)

# From :cite:`Iva99`, Lemma 4.4.1
DATA_GEOMETRY = {
   1: 48, 
   7: 759 * 8 * 2**7,
  22: binom(24,3) * 3 * 2**12,
   4: binom(24,4) * 2**4,
  10: 759 * binom(16,2) * 2**9,
  13: 2576 * 12 * 2**12,
  19: binom(24,5) * 2**12,
  16: 759 * 16 * 2**11,
}  
 
NUM_LEECH_TYPE2 = 398034000 
assert sum(DATA_GEOMETRY.values()) ==  NUM_LEECH_TYPE2
I_NUMV3 = 3.0**(-24)

DICT_P = defaultdict(int)
P = defaultdict(float)
for w, num in DATA_GEOMETRY.items():
    p = num * I_NUMV3
    assert 0 <= p < 1
    DICT_P[w] = 1 if p < 1.0/4000 else w
    P[DICT_P[w]] += p
P[0] = 1.0 - NUM_LEECH_TYPE2 * I_NUMV3
DICT_P[0] = 0

P_MIN = min([x for x in P.values() if x > 0])

RANDMOD3 = [0, 1, 0x1000000]
def rand_v3():
    randomlist = choices(RANDMOD3, k=24)
    return sum((x << i for i, x in enumerate(randomlist)))
    
def rand_v3_dict(n):
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
def test_chisq_type4(n = 50000, verbose = 0):
    p_min = 0.01
    d = rand_v3_dict(n)  
    print("Check distribution of type-4 vectors mod 3") 
    for i in range(4):
        chisq, p =  chisquare_v3(d, P)
        if verbose or i or p < p_min:
            print("Chisq = %.3f, p = %.4f" % (chisq, p))
        if p >= p_min: return
    raise ValueError("Chisquare test failed") 

#*************************************************************************
#** Self test from C file
#************************************************************************/


def one_selftest_leech2(data):
    start, n = data
    a = np.zeros(0x50, dtype = np.uint32)
    result =  gen_leech2_type_selftest(start, n, a)
    return result, a


def gen_selftest_inputs(n):
    assert 0x1000000 % n == 0
    q = 0x1000000 // n
    for  i in range(n):
        yield i*q, q



# From :cite:`Iva99`, Lemmas 4.4.1 and 4.6.1
TYPE_LENGTHS = {               # Name in :cite:`Iva99`
 0x00: 1,
 0x20: binom(24,2) * 2,        # \Lambda_2^4
 0x21: 24 * 2**11,             # \Lambda_2^3
 0x22: 759 * 2**6,             # \Lambda_2^2
 0x31: 24 * 2**11,             # \Lambda_3^5
 0x33: binom(24,3) * 2**11,    # \Lambda_3^3
 0x34: 759 * 16 * 2**7,        # \Lambda_3^4
 0x36: 2576 * 2**10,           # \Lambda_3^2
 0x40: 2 * 1771,               # \bar{\Lambda}_4^{4a}
 0x42: 759 * 2**6,             # \bar{\Lambda}_4^{6} 
 0x43: binom(24,3) * 2**11,    # \bar{\Lambda}_4^{5}
 0x44: 15 * 759 * 2**7,        # \bar{\Lambda}_4^{4b}
 0x46: 1288 * 2**11,           # \bar{\Lambda}_4^{4c}
 0x48: 1                       # \bar{\Lambda}_4^{8}
}  

@pytest.mark.gen_xi
@pytest.mark.slow
def test_leech2_self(verbose = 0):
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







