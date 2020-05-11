from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import pytest

import types
import sys
import re
import os
from operator import __or__, __xor__, __and__
from functools import reduce


 
from mmgroup.bitfunctions import v2, lmap, lrange

from mmgroup.generate_c  import prepend_blanks
from mmgroup.dev.mat24.mat24aux import generate_golay24_decode


def lsbit24(x):
    return v2(x | 0x1000000)


from mmgroup.dev.mat24.mat24heptad import HeptadCompleter
from mmgroup import mat24
from mmgroup.dev.mat24.mat24_ref import Mat24

def random_umbral_heptad(gc, as_heptad_at_pos = None):
    """Return random list p with p[0],...,p[5], p[8] an umbral heptad.

    That list hat length 24 (with unused entries set to None) and 
    it can be completed to a unique permutation in Mat24.

    So p[0],...,p[5]  are contained in an octad. The special element
    p[8] is not contained in that octad.

    if 0 <= as_heptad_at_pos < 7 is given then the heptad is returned
    as a list of 7 elements with the special element at that position.

    Parameter gc must be an instance of class Mat24.     
    """
    import random
    p = []
    while len(p) < 5:
        x = random.randint(0,23)
        if not x in p:
            p.append(x)
    as_cocode = gc.vect_to_vintern(sum(1 << x for x in p))
    syn = gc.syndrome_table[ as_cocode & 0x7ff ]
    syn = [ (syn >> i) & 31 for i in (0, 5, 10) ]
    octad, p8 = p + syn, None
    while p8 is None:
        x = random.randint(0,23) 
        if not x in octad:
            p8 = x
    p.append(syn[random.randint(0,2)])
    if not as_heptad_at_pos is None:
        return p[:as_heptad_at_pos] + [p8] + p[as_heptad_at_pos:]         
    return p + [None, None, p8] + [None]*15
     
   


def get_testdata(gc):
    hc = gc.heptad_completer
    import random
    prod = lrange(24)
    yield prod, True
    for i in range(100):
        p = random_umbral_heptad(gc)
        yield p, True
        hc.compute(p) 
        prod =  gc.mul_perm(prod, p)
        #print( prod )
        yield p, True

    # test some error cases
    for i in range(100):
        p = random_umbral_heptad(gc)
        hc.compute(p)
        p6 = p[6]
        j = random.randint(0,5)
        p[j], p[8] = p[8], p[j]
        yield p,  False
        p[j], p[8] = p[8], p[j]
        yield p,  True
        ok_j = p[j]
        p[j] = p[8]
        yield p,  False
        p[j] = random.randint(24,32)
        yield p,  False
        p[j] = ok_j
        yield p,  True
        p[8] = p6
        yield p,  False
        p[8] = p[j]
        yield p, False
        p[8] = random.randint(24,32)
        yield p,  False
    for i in range(6):
         for j in range(i):
             p = random_umbral_heptad(gc)
             hc.compute(p)
             p[i] = p[j]
             yield p, False     
              
        
  
def heptad_testdata(gc):
    hc = gc.heptad_completer
    import random    
    pos = random.randint(0,6)    
    for i in range(50):
        h1 = random_umbral_heptad(gc, pos)
        h2 = random_umbral_heptad(gc, pos)
        yield h1, h2


@pytest.mark.mat24
@pytest.mark.parametrize("gc, gc_ref", [(mat24, Mat24)])
def test_heptad_completer(gc, gc_ref):
    hc = gc_ref.heptad_completer
    for p, ok in get_testdata(gc_ref):
        p_in = p[:6] + [None]*2 + p[8:9] + [None]*15
        err =  gc.perm_complete_heptad(p_in)
        if ok:
            assert  err == 0, (p_in, hex(err))
            if not None in p:
                assert list(p) == list(p_in), (p, p_in, hex(err))
                i = hc.perm_to_int(p)
                assert 0 <= i < 244823040
                p1 = hc.int_to_perm(i)
                assert list(p) == list(p1), (i, p, p1)
        else:
            assert err != 0

    for h1, h2 in heptad_testdata(gc_ref):
        p = gc.perm_from_heptads(h1, h2)
        for i in range(7):
            assert p[h1[i]] == h2[i]

    perm0 = gc.perm_to_m24num(lrange(24))
    assert perm0 == 0, perm0

    print("Heptad completer test passed")











