from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import pytest

import types
import sys
import re
import os
from operator import __or__, __xor__, __and__
from functools import reduce
import random


 
from mmgroup.bitfunctions import v2, lmap, lrange

from mmgroup.generate_c  import prepend_blanks
from mmgroup.dev.mat24.mat24aux import generate_golay24_decode


def lsbit24(x):
    return v2(x | 0x1000000)


from mmgroup.dev.mat24.mat24heptad import HeptadCompleter
from mmgroup import mat24
from mmgroup.dev.mat24.mat24_ref import Mat24


#########################################################################
# Test completion of heptads and octads
#########################################################################


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
        hc.complete_heptad(p) 
        prod =  gc.mul_perm(prod, p)
        #print( prod )
        yield p, True

    # test some error cases
    for i in range(100):
        p = random_umbral_heptad(gc)
        hc.complete_heptad(p)
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
            hc.complete_heptad(p)
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
        if ok:
            err =  gc.perm_complete_heptad(p_in)
            #assert  err == 0, (p_in, hex(err))
            if not None in p:
                assert list(p) == list(p_in), (p, p_in, hex(err))
                i = hc.perm_to_int(p)
                assert 0 <= i < 244823040
                p1 = hc.int_to_perm(i)
                assert list(p) == list(p1), (i, p, p1)
        else:
            with pytest.raises(ValueError):
                gc.perm_complete_heptad(p_in)
        if ok:
            hexad = p[:6] + [None]*2
            hexad_ref = hexad[:]
            gc.perm_complete_octad(hexad)
            assert hexad == p_in[:8]
            gc_ref.perm_complete_octad(hexad_ref)
            assert hexad_ref == p_in[:8]

    for h1, h2 in heptad_testdata(gc_ref):
        p = gc.perm_from_heptads(h1, h2)
        for i in range(7):
            assert p[h1[i]] == h2[i]
    perm0 = gc.perm_to_m24num(lrange(24))
    assert perm0 == 0, perm0

    print("Heptad completer test passed")



#########################################################################
# Test completion of heptads and dodecads
#########################################################################

def random_dodecad():
    while True:
        gc = random.randint(0,0xfff)
        if mat24.gcode_weight(gc) == 3:
            v = mat24.gcode_to_vect(gc)
            w, l = mat24.vect_to_bit_list(v)
            assert w == 12
            l = l[:12]
            random.shuffle(l)
            return l


def get_test_dodecads():
    for i in range(200):
        yield random_dodecad(), random_dodecad()


@pytest.mark.mat24
def test_perm_from_dodecad():
    for ntest, (d1, d2) in enumerate(get_test_dodecads()):
        #print(d1, d2)
        perm = mat24.perm_from_dodecads(d1[:9], d2[:9])
        for i in range(5):
             assert perm[d1[i]] == d2[i]
        img = set([perm[x] for x in d1])
        assert img == set(d2)

    

#########################################################################
# Test function  mat24_perm_from_map()
#########################################################################



def rand_perm():
    """Return random permutation in ``Mat24`` as a list"""
    num = random.randint(0, mat24.MAT24_ORDER-1)
    return mat24.m24num_to_perm(num)


def map_vect(v):
    """Apply radom permutaiton in in ``Mat24`` of vector

    Here ``v`` is a vector of integers ``0 <= v[i] < 24``.
    The function generates a random permutation ``p`` in
    the group ``Mat24`` and returns the list
    
    [p[0], p[1],...,p[22], p[23]].
    """
    pi = rand_perm()
    return [pi[x] for x in v]


def perm_from_map_testdata():
    """Yield test data for testing function perm_from_map()

    The function yields quadruples (h1, h2, ref_res, ref_p).

    The function perm_from_map() is expected to find a
    permutation p that maps the entries in the list h1 to
    the entries in the list h2, provided that ``ref_res``
    is greater than zero.

    Also, function perm_from_map() should return the value
    ``ref_res``. If ``ref_p`` is not None then the permutation
    returnd by function perm_from_map() should be equal
    to ``ref_p``.    
    """
    Id = list(range(24))
    yield list(range(5)), list(range(5)), 3, Id
    yield [], [], 3, Id
    yield list(range(9)), list(range(9)), 1, Id
    yield [0,1,2,3,4,5,6], [0,1,2,3,4,5,8], 0, None

    PERMS = [
         ([0,1,2,3,4,9], 2),
         ([0,1,2,3,4,5,9], 1),
         ([0,1,2,3,4,5,6,9], 1),
    ]
    for i in range(1, 9):
        PERMS.append( (range(i), 3) )
    for perm, res in PERMS:
        v = map_vect(perm)
        yield v, v, res, Id
        for i in range(3):
            v1, v2 = map_vect(perm), map_vect(perm)
            yield v1, v2, res, None
            
    for n in range(9, 25):
        L = random.sample(range(24),n)
        yield L, L, 1, Id
        for j in range(3):
            p = rand_perm()
            ind = random.sample(range(24), n)
            perm = [p[i] for i in ind]
            yield ind, perm, 1, p
    

def one_test_perm_from_map(h1, h2, ref_res, ref_p,verbose = 1):
    if verbose:
        print("h1 =", h1)
        print("h2 =", h2)
        print("Expected return value:", ref_res)
    res, perm = mat24.perm_from_map(h1, h2)
    ok, ok_perm = res == ref_res, True
    if ref_p:
        ok_perm = perm == ref_p
    elif ref_res > 0:
        for i, x in enumerate(h1):
            ok_perm = ok_perm and perm[x] == h2[i]
    ok = ok and ok_perm
    if not verbose and not ok:
        print("h1 =", h1)
        print("h2 =", h2)
        print("Expected return value:", ref_res)
    if verbose or not ok:
        print("Obtained return value:", res)
        if res > 0:
            print("p =", perm)
            if ref_p and ref_p != perm:
                print("Expected:\n   ", ref_p)
        if not ok_perm:
            print("Permutation p is bad")
        if res != ref_res:
            print("Return value is bad")
    if not ok:
        err = "Test of function mat24.perm_from_map failed"
        raise ValueError(err)

    if ref_p or len(h1) <= 1:
        return

    # If no expected permutation is given, shuffle h1 and h2
    # in the same way and test once more with shuffled h1, h2.
    lh = len(h1)
    ind = list(range(lh)) if lh > 2 else [1,0]
    if (lh > 2):random.shuffle(ind)
    h1a = [h1[i] for i in ind]
    h2a = [h2[i] for i in ind]
    res1, perm1 = mat24.perm_from_map(h1a, h2a)
    assert res1 == res, (res1, res)
    assert perm1 == perm, (perm1, perm)


@pytest.mark.mat24
def test_perm_from_map(verbose = 0):
    for n, (h1, h2, ref_res, ref_p) in enumerate(perm_from_map_testdata()):
        if verbose: print("\nTest", n+1)
        one_test_perm_from_map(h1, h2, ref_res, ref_p, verbose)   
    