


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
from mmgroup.dev.generators.gen_leech_reduce_mod3 import reduce_leech_mod3
from mmgroup.dev.generators.gen_leech_reduce_mod3 import bitvector_syndromes
from mmgroup.dev.generators.gen_leech_reduce_mod3 import reduce_leech_mod3_testdata
from mmgroup.dev.generators.gen_leech_reduce_mod3 import test_reduce_leech_mod3


from mmgroup.generators import gen_leech3_find_tetrad_leech_mod3
from mmgroup.generators import gen_leech3_reduce_leech_mod3
from mmgroup.generators import gen_leech3_reduce


#####################################################################
# Test function gen_leech3_find_tetrad_leech_mod3
#####################################################################
#####################################################################
# Testing function reduce_type_2



def reduce_leech_mod3_C(a):
    ERR = "Function gen_leech3_reduce_leech_mod3 failed with error %s"
    g = np.zeros(12, dtype = np.uint32)
    res = gen_leech3_reduce_leech_mod3(a, g)
    if res < 0:
        a_red = gen_leech3_reduce(a)
        print("Input was 0x%012x" % gen_leech3_reduce(a_red))
        supp = (a_red | (a_red >> 24)) & 0xffffff
        w_v, w_syn, syn_l = bitvector_syndromes(supp)
        print("Weight: %d, syn_weight: %d" % (w_v, w_syn))
        for w_add, syn in syn_l:
            w_gc = w_v + w_syn - 2 * w_add 
            print("syn: 0x%06x, sub weight %d, gc: %d:" % 
                       (syn, w_add, w_gc))
        print(ERR % hex(abs(res)))
        raise ValueError(ERR % hex(abs(res)))
    return g[: res >> 48], res & 0xffffffffffff

@pytest.mark.mmm
@pytest.mark.gen_xi
def test_gen_leech3_find_tetrad_leech_mod3(verbose = 3):
    """Test function ``reduce_type2`` """
    print("\nTesting function gen_leech3_find_tetrad_leech_mod3")
    MSG = "Test %d, a = 0x%012x, tetrad expected: 0x%06x, obtained: 0x%06x" 
    ERR = "Test %d, a = 0x%012x, tetrad expected: 0x%06x, status: %d" 
    for n, a in enumerate(reduce_leech_mod3_testdata(10000)):
        a_red = gen_leech3_reduce(a)
        tet_ref = find_tetrad_leech_mod3(a)
        tet = gen_leech3_find_tetrad_leech_mod3(a)
        ok = tet == tet_ref
        g_ref, a1_ref = reduce_leech_mod3(a, verbose)
        g, a1 = reduce_leech_mod3_C(a)
        # Activate and test the following statment!!!!!!
        # ok = ok and list(g) == list(g_ref) and a1_ref == a1

        if verbose or not ok:           
            msg = MSG if tet >= 0 else ERR
            print(msg % (n+1, a_red, tet_ref, tet))
            a_red = gen_leech3_reduce(a)
            supp = (a_red | (a_red >> 24)) & 0xffffff
            neg = supp & (a_red >> 24)
            print("Support: 0x%06x, neg: 0x%06x" % (supp, neg))
            if not ok:
                w_v, w_syn, syn_l = bitvector_syndromes(a_red)
                print("Weight: %d, syn_weight: %d" % (w_v, w_syn))
                for w_add, syn in syn_l:
                    w_gc = w_v + w_syn - 2 * w_add 
                    print("syn: 0x%06x, sub weight %d, gc: %d:" % 
                       (syn, w_add, w_gc))
                msg = "Reduced word expected: 0x%012x, obtained: 0x%012x"
                print(msg % (a1_ref, a1))
                err = "Function en_leech3_find_tetrad_leech_mod3 failed"
                print("g expected", [hex(x) for x in g_ref])
                print("g obtained", [hex(x) for x in g])
                raise ValueError(err)
                pass

