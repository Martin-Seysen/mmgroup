"""Test C function dealing with Leech lattice vectors mod 3 of type 3

In this script we test functions dealing with vectors of type 3 in the 
Leech lattice modulo 3. These functions are implemented in file 
gen_leech3.c and available in the extension mmgroup.generators.

We use the terminology defined in
the document *The C interface of the mmgroup project*, 
section *Description of the mmgroup.generators extension*.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral
from collections import defaultdict

import numpy as np
import pytest

from mmgroup import MM0, Cocode, XLeech2
from mmgroup.generators import gen_leech2_op_word_leech2
from mmgroup.generators import gen_leech2_op_word_matrix24
from mmgroup.generators import gen_leech2_op_mul_matrix24



#####################################################################
# Creating test elements of G_x0 and of leech lattice mod 2
#####################################################################

def Co_1_samples():
    testdata = [
      (0x800000,   [], []),
      (randint(0, 0xfffff),  [('p', 1000)],  [('y', 0x123)]),
      (randint(0, 0xfffff),  [('l', 1)],  [('l', 2)]),
    ]
    for x, d1, d2  in testdata:
        yield x,   MM0(d1).mmdata, MM0(d2).mmdata
    for i in range(100):
        yield (randint(0, 0xffffff), MM0('r', 'G_x0').mmdata, 
            MM0('r', 'G_x0').mmdata)







#####################################################################
# Test mapping of type 3 Leech vectors modulo 3 to vectors modulo 2
#####################################################################



preimages = [0x800000] + [Cocode([0,i]).ord for i in range(1,24)]



@pytest.mark.gen_xi
def test_op_word_matrix24(verbose = 0):
    """Test the C function ``gen_leech2_map_std_subframe``"""
    print("Testing  C function gen_leech2_op_word_matrix24")
    a1 = np.zeros(24, dtype = np.uint32)
    a2 = np.zeros(24, dtype = np.uint32)
    y = np.zeros(1, dtype = np.uint32)
    for x, g1, g2 in Co_1_samples():
        if verbose:
           print("Testing %s * %s * %s**-1" % (XLeech2(x),
               MM0('a', g1), MM0('a', g2)))
        assert gen_leech2_op_word_matrix24(g1, len(g1), 0, a1) >= 0
        assert gen_leech2_op_word_matrix24(g2, len(g2), 1, a2) >= 0
        #print(a1, a2)
        x1 =  gen_leech2_op_word_leech2(x, g1, len(g1), 0)
        x2 =  gen_leech2_op_word_leech2(x1, g2, len(g2), 1)
        gen_leech2_op_mul_matrix24(a1, 24, a2)
        #print(a1)
        y[0] = x
        if verbose:
           print("Result: %s" % hex(x2))
        gen_leech2_op_mul_matrix24(y, 1, a1)
        assert y[0] == x2, (hex(x), hex(x2), hex(y[0]))
    print("passed")
                 



