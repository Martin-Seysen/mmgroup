from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


from numbers import Integral
import numpy as np
from random import randint
import math

import pytest

from mmgroup.mm_crt_space import MMSpaceCRT
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
    yield([(tag,)])

def vector(tags, max_norm):
    res = []
    sc = abs_rand_int(max_norm)
    for tag in tags:
        res.append((sc, tag, "r"))
    yield res
    

def crt_test_tuples(k):
    max_norm = 1.234 * 10**13 * 4.0 ** (-k)
    for tag in "BCTZXYDE":
        yield from units(tag, math.sqrt(max_norm))
    yield from units("A", math.sqrt(max_norm/2))
    yield from units("I", math.sqrt(max_norm/8))
    for tag1 in "ABCTXZY":
         for tag2 in "ABCTXZY":
             if tag1 < tag2:
                 yield from vector(tag1+tag2, math.sqrt(max_norm/3))
             
def crt_testdata(k):
    group = standard_mm_group
    yield [("A",0,0)],  group(("d", 0x801))
    yield [("A",1,0)],  group(("t", 1))
    for g_tags in ["tlxypd" * 2, "tlxypd" * 4 ]:
         for v in crt_test_tuples(k): 
             g = group.sample(g_tags, None, "n") 
             yield v, g  



def g_shift(g):
    shift = 0
    for g0 in g.data:
        if g0 & 0x70000000 >= 0x50000000:
            shift += 3
    return shift

@pytest.mark.mm
def test_random_io(verbose = 0):
    n = 1
    for k in (17, 20):
        space =  MMSpaceCRT(k)
        for v, g in crt_testdata(k):
             if verbose:
                 print("\nTest", n)
                 print("v =", v)
                 print("g =", g)
             v = space(*v)
             n1 = v.inorm() 
             v2_1 = v.v2()
             v *= g     
             n2 = v.inorm() 
             assert n1 == n2, (hex(n1), hex(n2))
             for p in (7,31, 127, 255):
                 vp =  (v % p) 
                 for tag in "A": # "ABCTXZY":
                     assert np.all(v[tag] % p == vp[tag])
             v2_2 = v.v2()
             assert v2_2 >= v2_1 - g_shift(g)
             n += 1


