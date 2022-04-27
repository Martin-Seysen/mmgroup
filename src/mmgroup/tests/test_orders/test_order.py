from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import numpy as np
from numbers import Integral
from random import randint, sample


import pytest

from mmgroup import  MMV, MM0

ref_space = MMV(3)
ref_group = MM0

from mmgroup.structures import mm_order
#from mmgroup.structures.mm_order import reduce_mm

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
   





