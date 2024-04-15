from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import numpy as np
from numbers import Integral
from random import randint, sample
import math

import pytest

from mmgroup import  MMV, MM0

ref_space = MMV(3)
ref_group = MM0

def import_all():
    global mm_order
    from mmgroup.structures import mm_order

########################################################################
# Computing the order of a group element using a random vector
########################################################################

def ref_order_random_v(g):
    g1 = ref_group(g)
    v = ref_space("R")
    w = v.copy()
    for i in range(1, 120):
        w = w * g1
        if w == v:
            return i
    err = "Error in computation reference order of group element"
    raise ValueError(err)
 
 
def gcd(a, b):
    a, b = max(a,b), min(a,b)
    while b > 0:
        a, b = b, a % b
    return a


def ref_order(g):
    """Return order of element ``g`` of the monster

    The result is based on calculations with random vectors in the 
    representation modulo 3. It is correct with very high probability.
    """
    order = 1
    for k in range(3):
        order1 = ref_order_random_v(g)
        order = order * order1 // gcd(order, order1)
        if order > 59:
            assert order < 120
            return order
    assert order < 120
    return order


########################################################################
# Test computation of the order of a group element
########################################################################

def order_testcases(group):
    yield group(), 1
    for x,y in zip(range(0, 0x2000, 0x800), [0, 0x1800, 0x800, 0x1000]):
        yield  group([('x', x), ('y', y)]), 1
    yield group('t', 1), 3
    tags_list = [ 
       ("txpyd", 10), ("txpyld"*5, 1), 
       ("txpyld"*3, 3),  ("lxypdt", 5), 
       ("xypd", 3),
    ]
    for i in range(2):
        for tags, ntests in tags_list:
            for j in range(ntests):
                g = group([(t,"n") for t in tags])
                order =  ref_order(g)
                yield g, order

   
@pytest.mark.orders
def test_order(verbose = 0):
    import_all()
    group = MM0
    for n, (g, ref_order) in enumerate(order_testcases(group)):
        order = g.order()
        ok = order == ref_order
        if verbose:
            print("Test", n + 1)
        if verbose or not ok:
            print("g =", g)
            print("order =", ref_order)
        if not ok:
            print("computed order =",order)
            err = "Error in computation order of group element"
            raise ValueError(err)
        ho, g2 = g.copy().half_order()
        assert ho == order
        if order & 1:
            assert g2 == None
        else:
            assert g2 == g**(order // 2)
    print("Test passed")
   
   


########################################################################
# Test computation of equality of group elements
########################################################################

def ref_equal_random(g1, g2):
    """ ``g1 == g2`` for elements ``g1, g2`` of the monster

    The result is based on calculations with random vectors in the 
    representation modulo 3. It is correct with very high probability.

    Elements of different instances of class MMGroup are 
    considered as unequal.
    """
    if g1.group != g2.group:
        return False
    g1 = ref_group(g1).reduce()
    g2 = ref_group(g2).reduce()
    g = g1 * g2**(-1)
    g.reduce()
    for i in range(3):
        v = ref_space('R')
        w = v.copy()
        w = w * g
        if w != v:
            return False
    return True
 





def equality_testcases(group):
    yield group(), group(), 1
    tags_list = [ 
       ("txpyd", 3), 
       ("txpyld"*2, 2), 
       ("xypdtl", 2), ("xypd", 3),
    ]
    for i in range(5):
        for tags, ntests in tags_list:
            for j in range(ntests):
                g = group([(t,"n") for t in tags])
                order =  g.order()
                assert g**order == group()
                prefix = group('r', randint(5,10))
                suffix = group('r', randint(5,10))
                prefix = group(1)
                suffix =  group(1)
                half = order//2
                g1 = prefix * g**(-half) * suffix
                g2 = prefix * g**(order - half) * suffix
                yield g1, g2, True
    for i in range(1):
        for tags, ntests in tags_list:
            for j in range(ntests):
                g1 = group([(t,"n") for t in tags])
                g2 = group([(t,"n") for t in tags])
                yield g1, g2, ref_equal_random(g1, g2)


@pytest.mark.orders
def test_equality(verbose = 0):
    import_all()
    group = MM0
    for n, (g1, g2, ref_equal) in enumerate(equality_testcases(group)):
        equal = g1 == g2
        ok = equal == ref_equal
        if verbose:
            print("Test", n + 1)
        if verbose or not ok:
            print("g1 =", g1)
            print("g2 =", g2)
            print("equal" if ref_equal else "unequal")
        if not ok:
            err = "Equality test of group elements failed"
            raise ValueError(err)
    print("Test passed")
   

################################################################################
# Testing method chi_powers() of class MM
################################################################################


MMGROUP = MM0                   # We compute in this base class of class MM
STD_2A = MMGROUP('d', [2,3])    # The standard 2A involution
STD_2B = MMGROUP('x', 0x1000)   # The standard 2B involution


def chi_power_samples():
    """Samples for testing method chi_powers() of class MMGROUP"""
    g_samples = [
        # Erroneous sample (in old version) found by Gerald Hoehn
        "M<y_754h*x_11fah*d_4c1h*p_167735899*l_1*p_1499520*l_1*p_33398856*l_2*t_2*l_1*p_24000*l_2*p_32018160*l_1*t_2*l_2*p_1394880*l_1*p_10668720*l_1*p_520320*t_2*l_2*p_3569280*l_1*t_1*l_2*p_2386560*l_2*p_21338083*t_2*l_1*p_1499520*l_2*p_96465200*t_1*l_2*p_48829440*l_2*p_85198272*t_1*l_2*p_2787840*l_2*p_64736>",
        # Basis of following sample in class 4B, see file involution_samples.py
        MMGROUP("M<y_2e0h*x_1020h*d_0fh>") ** MMGROUP('r', 2),
        # Conjugates of standard involutions
        STD_2A ** MMGROUP('r', 2),
        STD_2B ** MMGROUP('r', 2),
    ]
    for g in g_samples:
        yield MMGROUP(g)
    for i in range(3):
        yield MMGROUP('r', 'G_x0')
"""
CHAR_M is the set of possible values of the character of the Monster.
Values of this character have been computed with the following GAP program:

LoadPackage( "AtlasRep", false );
tbl:= CharacterTable( "M" );;
Irr(tbl)[2];
"""
CHAR_M = set([
  196883, 4371, 275, 782, 53, -1, 275, 51, 19, -13, 133, 8, 78, 77, 14, -3, 5, -1,
  50, 1, 35, 11, -1, -5, 3, -1, 26, -1, 21, 5, -4, 20, 0, 16, 14, 5, 6, -1, -2, 5, -3, 13, 1, -1, 11, -2, 10, 2, 9,
  7, -2, 8, -1, 3, -1, 7, 6, -3, 6, 2, -1, 5, 5, 5, 1, 0, -3, 2, 4, 5, -2, -1, 4, 4, 0, 3, 3, 2, 2, -1, -2, -1, -1,
  -1, 1, 3, -1, 3, 3, 2, 2, 2, 2, 2, -2, 1, 2, 2, 3, -1, 2, -1, 2, 0, 2, 2, 1, 1, -2, 1, 2, 0, 1, 2, -1, 0, 1,
  1, 2, -1, 1, 1, -1, 1, 0, 0, 1, 1, 0, -1, 0, 0, 0, 1, -1, -1, 1, 1, 0, 0, 0, 1, 0, -1, 0, 0, 1, 0, -1, -1, -1, 0,
  0, 1, -1, 0, -2, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, -1, -1, -2, -1, -1, -1, 0, 0, -1, -1, -1, -1, 0,
  0, 0, 0, -1, -1, 0, -1, -1, -1 
])

# Add None to the set CHAR_M, for the case that method chi_powers() fails
CHAR_M.add(None)

def set_of_divisors(order):
    divisors = set([1, order])
    for i in range(2, int(math.sqrt(order) + 1.01)):
        if order % i == 0:
            divisors |=  set([i, order // i])
    return divisors
    
@pytest.mark.orders
def test_chi_powers():
    """Test method g.chi_powers() of class MMGROUP"""
    import_all()
    print("")
    MAX_E = None
    found_2A = found_2B = False
    for g in chi_power_samples():
        assert isinstance(g, MMGROUP)
        order, chi, h = g.chi_powers(max_e=MAX_E)
        assert order == g.order()
        assert chi.keys() == set_of_divisors(order)     
        for chi_e in chi.values():
            assert chi_e in CHAR_M
        assert isinstance(h, MMGROUP)
        if order & 1 == 0:
            chi_sqrt_1 = chi[order//2]
            assert chi_sqrt_1 in [None, 275, 4371]
            g1 = g**h
            if chi_sqrt_1 == 275:  # g is in class 2B
                found_2B = True
                assert g1.in_G_x0()
                assert 1 in chi and 2 in chi
                for e, chi_e in chi.items():
                    if MAX_E and e <= MAX_E:
                        assert chi_e == (g1**e).chi_G_x0()[0]
            if chi_sqrt_1 == 4371:  # g is in class 2A
                found_2A = True
                assert  STD_2A ** g1 == STD_2A

    assert found_2A and found_2B






