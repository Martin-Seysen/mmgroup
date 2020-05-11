from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import os
import collections
import re
import warnings
from numbers import Integral
import numpy as np

import pytest

from mmgroup.tests.groups.mgroup_n import MGroupN
from mmgroup.tests.groups.mgroup_n import n_rules, n_inv


###########################################################################
# Test functions
###########################################################################

MIN_VERBOSE = 2




RAND_TAGS = "dpxyt" * 10
N_TESTS = 200



@pytest.mark.auto_group
def test_mgroupn(n_tests = N_TESTS, rand_tags = RAND_TAGS, verbose = 0):
    global n_rules, n_inv
    n_rules =  n_inv = 0
    g = MGroupN(0)
    x1 = g.atom('x', 0x1234)
    x2 = g.atom('x', 0x1fe7)
    y1 = g.atom('y', 0x1fe1)
    y2 = g.atom('y', 0x1e31)
    p = g(('d', 0x3b5), ('p', 123423458))
    print("x1*x2",  x1, x2,  (x1 * x2)**(-3) )
    print( g.sample(rand_tags).reduce() )
    for i in range(n_tests):
        #w1, w2, w3 = g(rand_tags), g(rand_tags), g(rand_tags)
        w1, w2 = g.sample(rand_tags), g.sample(rand_tags)
        w3 = g.sample(rand_tags) 
        assert g(str(w1)) == w1, (g(str(w1)), w1)
        assert g(*w1.as_tuples()) == w1
        l_assoc = ((w1 * w2) * w3).reduce()
        r_assoc = w1 * (w2 * w3)
        ok = l_assoc == r_assoc
        if verbose or not ok or i < MIN_VERBOSE:
            print("\nTest associative law")
            print("w1 =", w1)			
            print("w2 =", w2)			
            print("w3 =", w3)			
            print("(w1*w2)*w3 =", l_assoc )
            if not ok:			
                print("w1*(w2*w3= =", w1 * r_assoc )			
                raise ValueError("Test of accociative law Failed")
        a = (w1 * w2 * w3)**(-1) 
        b =  w3**(-1) * w2**(-1) * w1**(-1)
        ok = a == b
        if verbose or not ok or i < MIN_VERBOSE:
            print("\nTest Inversion of w1 * w2 * w3")
            print("w1 =", w1)			
            print("w2 =", w2)
            print("w3 =", w3)
            print("a = (w1 * w2 * w3)**(-1)")			
            print("a =", a)			
            print("b = w3**(-1) * w2**(-1) * w1**(-1)")			
            print("b =", b)
            if not ok:			
                print("a * b**(-1)) =", a * b**(-1)) 
                print("lengths of a and b",  len(a), len(b))
                print("atom comparison", [ x == y for x, y in zip(a, b)] )
                raise ValueError("Inversion Failed")
    print("%d rules, %d inversions" % (n_rules, n_inv))





#if __name__ == "__main__":
#    test_mgroupn()


