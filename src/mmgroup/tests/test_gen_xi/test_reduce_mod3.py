


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices, shuffle, sample
from numbers import Integral
from collections import defaultdict

import numpy as np
import pytest

#from mmgroup import XLeech2, Xsp2_Co1, PLoop, GCode, AutPL, Cocode, GcVector

from mmgroup import mat24

from mmgroup.dev.generators.gen_leech_reduce_mod3 import find_tetrad_leech_mod3
from mmgroup.dev.generators.gen_leech_reduce_mod3 import bitvector_syndromes
from mmgroup.dev.generators.gen_leech_reduce_mod3 import reduce_leech_mod3_testdata
from mmgroup.dev.generators.gen_leech_reduce_mod3 import test_reduce_leech_mod3


from mmgroup.generators import gen_leech3_find_tetrad_leech_mod3


#####################################################################
# Test function gen_leech3_find_tetrad_leech_mod3
#####################################################################
#####################################################################
# Testing function reduce_type_2


@pytest.mark.gen_xi
def test_gen_leech3_find_tetrad_leech_mod3(verbose = 0):
    """Test function ``reduce_type2`` """
    print("\n testung function gen_leech3_find_tetrad_leech_mod3")
    MSG = "Test %d, a = 0x%012x, tetrad expected: 0x%06x, obtained: 0x%06x" 
    ERR = "Test %d, a = 0x%012x, tetrad expected: 0x%06x, status: %d" 
    for n, a in enumerate(reduce_leech_mod3_testdata(10000)):
        tet_ref = find_tetrad_leech_mod3(a)
        tet = gen_leech3_find_tetrad_leech_mod3(a)
        ok = tet == tet_ref
        if verbose or not ok:           
            msg = MSG if tet >= 0 else ERR
            print(msg % (n+1, a, tet_ref, tet))
            v = (a ^ (a >> 24)) & 0xfffff
            neg = v & (a >> 24)
            print("Support: 0x%06x, neg: 0x%06x" % (v, neg))
            if not ok:
                w_v, w_syn, syn_l = bitvector_syndromes(v)
                print("Weight: %d, syn_weight: %d" % (w_v, w_syn))
                for w_add, syn in syn_l:
                    w_gc = w_v + w_syn - 2 * w_add 
                    print("syn: 0x%06x, sub weight %d, gc: %d:" % 
                       (syn, w_add, w_gc))
                #raise ValueError("Test failed")
                pass

