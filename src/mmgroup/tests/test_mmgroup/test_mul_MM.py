
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import sample
import datetime
import time


import pytest

from mmgroup import MM, MMV

#####################################################################################
# Test fast reduction in monster group
#####################################################################################

MIN_LEN_ALWAYS_UNREDUCED = 80

def mul_testcases(ncases = 1):
    yield MM(1)
    for tag in "yxdptl":
         yield MM(tag, 'r')
    for complexity in range(8):
         for i in range(ncases):
              yield MM('r', complexity)

@pytest.mark.mmgroup 
def test_mul_MM(verbose = 0):
    print("Test multiplcation in monster group in class MM")
    TEXT = "Test %d: multiply random vector v (mod 3) with g1 and g2"
    V = MMV(3)
    n = 1
    for g1 in mul_testcases():
        v = V("R")
        for g2 in  mul_testcases():
            if verbose:
                print(TEXT % n)
                print(" g1 =", g1)
                print(" g2 =", g2)
                print(" Computing (v * g1) * g2")
            v1 = (v * g1) * g2
            v2 = v * (g1 * g2)
            if verbose:
                print(" Computing  v * (g1 * g2)")
            assert v1 == v2
            if verbose:
                print(" passed")
                n += 1

