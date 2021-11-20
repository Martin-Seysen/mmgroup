from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices #, shuffle, sample
from numbers import Integral
from collections import defaultdict

import numpy as np
import pytest

from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_op_word
from mmgroup import XLeech2, Xsp2_Co1, PLoop, GCode, AutPL, Cocode, GcVector
from mmgroup.generators import gen_leech2_start_type4
from mmgroup.generators import gen_leech2_start_type24


from mmgroup.tests.test_gen_xi.reduce import leech2_start_type4




OMEGA = 0x800000
STD_V2 = 0x200 

assert Cocode(STD_V2).syndrome() == GcVector(0xC)

#####################################################################
# Test function gen_leech2_start_type4
#####################################################################

def ref_leech2_start_type4(v):
    v &= 0xffffff
    t = gen_leech2_type(v)
    if (t & 0xf0) != 0x40:
        return min(-1, -(t >> 4))
    if v == OMEGA:
        return 0
    t2 = gen_leech2_type(v ^ STD_V2)
    return t2 if (t2 & 0xf0) == 0x20 else t

def type4_testdata(ntests = 1000):
    testdata = [
        0, 
        0x800,
        0xC,
    ]
    for t in testdata:
         yield t
    for type_ in range(2,4):
        for i in range(5):
            yield XLeech2('r', type_).ord 
    for i in range(24):
         yield (1 << i) + ((randint(1,10) << 24))
         yield (1 << i) ^ 0xffffff
    for i in range(ntests):
        yield XLeech2('r', 4).ord
    for i in range(100):
        v =  XLeech2('r', 2) * XLeech2(STD_V2)
        if v.type == 4:
            return v.ord

@pytest.mark.xsp2co1
def test_start_type4(ntests = 1000, verbose = False):
          
    for n, v in enumerate(type4_testdata(ntests)):
        if verbose:
            print("Test %d, v = %s, subtype = %s, subtpe(v2) = %s" % (
                n, hex(v), hex(gen_leech2_type(v)),
                hex(gen_leech2_type(v ^ STD_V2)) 
        ))
        t_py = leech2_start_type4(v)
        t_ref = ref_leech2_start_type4(v) 
        assert t_py == t_ref, (hex(t_py), hex(t_ref))
        t_c = gen_leech2_start_type4(v)
        assert t_c == t_ref, (hex(t_c), hex(t_ref))
    


#####################################################################
# Test function gen_leech2_start_type24
#####################################################################


def ref_leech2_start_type24(v):
    v &= 0xffffff
    t = gen_leech2_type(v)
    if (t & 0xf0) != 0x20:
        return -1
    if (v & 0x7fffff) == 0x200:
        return 0 if v & 0x800000 else -1
    t2 = gen_leech2_type(v ^ STD_V2)
    return t if (t2 & 0xf0) == 0x40 else -1


def type24_testdata(ntests):
    XLeech2_STD_V2 = XLeech2(STD_V2)
    testdata = [
        0, 
        0x800,
        0xC,
        0x200,
        0x800200,
    ]
    for t in testdata:
         yield t
    for type_ in range(2,5):
        for i in range(5):
            yield XLeech2('r', type_).ord 
    for i in range(24):
         yield (1 << i) + ((randint(1,10) << 24))
         yield (1 << i) ^ 0xffffff
    for i in range(ntests):
        yield XLeech2('r', 2).ord
    for i in range(200):
        v =  XLeech2('r', 2)
        if (v * XLeech2_STD_V2).type == 4:
            return v.ord



@pytest.mark.xsp2co1
def test_start_type24(ntests = 1000, verbose = 0):          
    for n, v in enumerate(type24_testdata(ntests)):
        if verbose:
            print("Test %d, v = %s, subtype = %s, subtpe(v2) = %s" % (
                n, hex(v & 0xffffff), hex(gen_leech2_type(v)),
                hex(gen_leech2_type(v ^ STD_V2)) 
        ))
        #t_py = leech2_start_type24(v)
        t_ref = ref_leech2_start_type24(v) 
        #assert t_py == t_ref, (hex(t_py), hex(t_ref))
        t_c = gen_leech2_start_type24(v)
        assert t_c == t_ref, (hex(t_c), hex(t_ref))
    

