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

from mmgroup import MM0
from mmgroup.mat24 import MAT24_ORDER, ploop_theta
from mmgroup.mat24 import bw24 as mat24_bw24
from mmgroup.generators import gen_leech3to2_type3
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.generators import gen_leech3_op_vector_atom
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_type_selftest


OMEGA0_V3  = 0x1






#####################################################################
# Creating test vectors
#####################################################################

STD_TYPE3_MOD3 = 0xfffffe000001
STD_TYPE3_MOD2 = 0x800800



def rand_xsp2co1_elem(n):
    s1, l = [('p','r'),('y','r')],  [('l','n')]
    return s1 + (l + s1) * n






def create_test_elements():
    group_data = [
      [],
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
# Test mapping of type 3 Leech vectors modulo 3 to vectors modulo 2
#####################################################################


def mul_v3(v3, g):
    result = gen_leech3_op_vector_word(v3, g._data, g.length)
    assert result & 0xffff000000000000 == 0, hex(result)
    return result
        
def mul_v2(v2, g):
    result = gen_leech2_op_word(v2, g._data, g.length)
    assert result & 0xfe000000 == 0, hex(result)
    return result

def v3_to_v2(v3):
    result = gen_leech3to2_type3(v3)
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
def test_type3(verbose = 0):
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
        #if  weights[w] <= 20:
        #     assert  v2 == py_gen_leech3to2_type3(v3)        
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
    assert set(weights.keys()) == set([9, 12, 21, 24])
    
    
    
#####################################################################
# Chisquare test of random type 3 Leech vectors modulo 3
#####################################################################



def binom(n, k):
    return int(scipy.special.binom(n, k) + 0.1)

# From :cite:`Iva99`, Lemma 4.4.1
DATA_GEOMETRY = {
  24: 24 * 2**12, 
   9: 759 * 16 * 2**8,
  21: binom(24,3) * 2**12,
  12: 2576 * 2**11,
}  
 
NUM_LEECH_TYPE3 = 2**24 - 2**12
assert sum(DATA_GEOMETRY.values()) ==  NUM_LEECH_TYPE3
I_NUMV3 = 3.0**(-24)

BLOCKSIZE = 1000000
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
    randomlist = choices(RANDMOD3, k=24)
    return sum((x << i for i, x in enumerate(randomlist)))
    
def rand_v3_dict(n = BLOCKSIZE):
    d = defaultdict(int)
    for i in range(n):
        v3 = rand_v3()
        v2 = gen_leech3to2_type3(v3) 
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
    print(expected_dict.keys())
    print(f_exp)
    print(f_obt)
    assert min(f_exp) > 9,  f_exp 
    chisq, p = chisquare(f_obt, f_exp = f_exp)
    return chisq, p




@pytest.mark.slow
@pytest.mark.very_slow
@pytest.mark.gen_xi
def test_chisq_type3(verbose = 0):
    p_min = 0.01
    print("Check distribution of type-3 vectors mod 3") 
    for i in range(4):
        d = rand_v3_dict()  
        chisq, p =  chisquare_v3(d, P)
        if verbose or i or p < p_min:
            print("Chisq = %.3f, p = %.4f" % (chisq, p))
        if p >= p_min: return
    raise ValueError("Chisquare test failed") 



