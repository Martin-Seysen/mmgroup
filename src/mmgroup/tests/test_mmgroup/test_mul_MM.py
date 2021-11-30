
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
def test_mul_MM():
    print("Test multiplcation in monster group in class MM")
    V = MMV(3)
    for g1 in mul_testcases():
        v = V("R")
        for g2 in  mul_testcases():
            v1 = (v * g1) * g2
            v2 = v * (g1 * g2)
            assert v1 == v2

