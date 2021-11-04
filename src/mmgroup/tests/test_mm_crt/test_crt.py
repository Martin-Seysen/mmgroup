from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


from numbers import Integral, Number
import numpy as np
from random import randint
import math

import pytest

from mmgroup.mm_crt_space import MMSpaceCRT, MMV_CRT, MMVectorCRT
from mmgroup.mm_space  import MMSpace, standard_mm_group



#########################################################################
### Test io an elementary operations
#########################################################################



def abs_rand_int(x):
    x = int(math.floor(x))
    return randint(1, x) * (-1)**randint(0,1)



def units(tag, max_norm):
    sc1 = abs_rand_int(max_norm)
    sc2 = abs_rand_int(max_norm)
    if tag in "ABCTXZYI":
         r1 = r2 = randint(0,23)
         while r1 == r2:
             r2 = randint(0,23)
         yield [(sc1, tag, r1, r2)]
         yield [(sc2, tag, "r", "r")]
    else:
         yield [(sc1, tag, randint(0,23))]
         yield [(sc2, tag, "r")]
    yield([(tag, "r")])

def vector(tags, max_norm):
    res = []
    sc = abs_rand_int(max_norm)
    for tag in tags:
        res.append((sc, tag, "r"))
    yield res
    

def crt_test_tuples(k):
    max_norm = MMVectorCRT.MAX_CRT_NORM * 4.0 ** (-k)
    #print (max_norm,  MMVectorCRT.MAX_CRT_NORM* 4.0 ** (-k))
    
    for tag in "BCTZXYDE":
        yield from units(tag, math.sqrt(max_norm/4))
    yield from units("A", math.sqrt(max_norm/2))
    yield from units("I", math.sqrt(max_norm/8))
    for tag1 in "ABCTXZY":
         for tag2 in "ABCTXZY":
             if tag1 < tag2:
                 yield from vector(tag1+tag2, math.sqrt(max_norm/3))
             
def crt_testdata(k):
    group = standard_mm_group
    yield [("A",0,0)],  group("d", 0x801)
    yield [("A",1,0)],  group("t", 1)
    for n in (2,3):
         for v in crt_test_tuples(k): 
             g = group('r', n) 
             yield v, g  


SH = 1 << 25
def equ_modp(a1, a2, p):
    if isinstance(a1, Number) and  isinstance(a2, Number):
         diff = (a1 * SH - a2 * SH) % p
         return diff == 0
    if a1.shape == a2.shape:
        b1, b2 = a1.astype(float), a2.astype(float)
        diff = (b1 * SH - b2 * SH) % p
        #print(p);print(a1); print(a2)
        return (diff == 0).all()
    raise TypeError("Arrays cannot be compared")
       


def g_shift(g):
    shift = 0
    for g0 in g.mmdata:
        if g0 & 0x70000000 >= 0x50000000:
            shift += 3
    return shift

@pytest.mark.slow
@pytest.mark.mm_op_crt
def test_random_io(verbose = 0):
    n = 1
    for k in [17, 20]:
        space =  MMV_CRT(k)
        for v, g in crt_testdata(k):
             if verbose:
                 print("\nTest", n, ", k=", k)
                 print("v =", v)
                 print("g =", g)
             v = space(v)
             n1 = v.inorm 
             v2_1 = v.v2
             v1 = v * g     
             n2 = v.inorm 
             assert n1 == n2, (hex(n1), hex(n2))
             for p in (15, 7, 31, 127, 255):
                 vp =  (v % p) * g
                 for tag in "BC":
                     i0 = i1 = randint(0,23)
                     while i1 == i0: i1 = randint(0, 23)
                     assert equ_modp(v1[tag,i0,i1], vp[tag,i0,i1], p)
                 for tag in "A":
                     assert equ_modp(v1[tag], vp[tag], p)
                 for tag in "TXZY":
                     i0, i1 = randint(0,758), randint(0,23)
                     assert equ_modp(v1[tag,i0,i1], vp[tag,i0,i1], p)
                 for j in range(2):
                     i = randint(0, 851)
                     assert equ_modp(v1["E", i], vp["E", i], p), (p,i) 
                 for j in range(5):
                     i = randint(852, 196883)
                     assert equ_modp(v1["E", i], vp["E", i], p), (p,i) 
             v2_2 = v1.v2
             assert v2_2 >= v2_1 - g_shift(g)
             i = randint(0, 23)
             assert v1["E",i] == v1["A", i, i]
             n += 1


