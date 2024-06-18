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

from mmgroup import MM0, Cocode
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_map_std_subframe




#####################################################################
# Creating test elements of G_x0
#####################################################################

def G_x0_samples():
    testdata = [
        [],
        [('p', 1000)],
        [('l', 1)],
    ]
    for d in testdata:
        yield MM0(d).mmdata
    for i in range(100):
        yield MM0('r', 'G_x0').mmdata






#####################################################################
# Test mapping of type 3 Leech vectors modulo 3 to vectors modulo 2
#####################################################################



preimages = [0x800000] + [Cocode([0,i]).ord for i in range(1,24)]



@pytest.mark.gen_xi
def test_maps_std_subframe(verbose = 0):
    """Test the C function ``gen_leech2_map_std_subframe``"""
    print("Testing  C function gen_leech2_map_std_subframe")
    a = np.zeros(24, dtype = np.uint32)
    a_ref = np.zeros(24, dtype = np.uint32)
    BAD_MASK = -0x2000000 & 0xffffffff
    for g in  G_x0_samples():
        if verbose:
           print("Input", MM0('a', g)) 
        for i, v in enumerate(preimages):
            a_ref[i] = gen_leech2_op_word(preimages[i], g, len(g))
            assert a_ref[i] & BAD_MASK == 0
        res = gen_leech2_map_std_subframe(g, len(g), a)
        assert res == len(g)
        ok = (a == a_ref).all()
        if verbose or not ok:
            print("Output obtained")
            print(a)
            print("Output expected")
            print(a_ref)
            if not ok:
               raise ValueError("Test failed")
    print("passed")
                 



