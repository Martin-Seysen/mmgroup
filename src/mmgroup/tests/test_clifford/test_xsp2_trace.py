from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from multiprocessing import Pool
import time

import numpy as np
import pytest

from mmgroup import Xsp2_Co1, PLoop, AutPL, Cocode, MM
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_is_type2
from mmgroup.generators import gen_leech2_count_type2
from mmgroup.generators import gen_leech2_op_word
from mmgroup.clifford12 import xsp2co1_leech2_count_type2
from mmgroup.clifford12 import xsp2co1_trace_98280

#####################################################################
# Create test matrices
#####################################################################



def rand_xsp2co1_elem():
    return Xsp2_Co1(*[(x,) for x in "dxpylpylpylpy"])

def xsp2xco1_xsp2(v):
    pl = PLoop(v >> 12)
    coc = Cocode(v)
    return Xsp2_Co1(("x", v >> 12), pl.theta(), coc)


v2_types = {0:0}


def xsp2xco1_v2type(vtype):
    global v2_types
    try:
        return xsp2xco1_xsp2(v2_types[vtype])
    except KeyError:
        assert vtype in [2,3,4] 
        for i in range(10000):
            v = randint(0, 0xffffff)
            if gen_leech2_type(v) >> 4 == vtype:
                v2_types[vtype] = xsp2xco1_xsp2(v) 
                return v2_types[vtype]
        raise ValueError("No Leech lattice vector of type", vtype)


@pytest.mark.xsp2co1a
def test_xsp2_count_table():
    use_table = 0
    table = [0] * 5
    for vtype in [0,2,3,4]:
        tr = np.zeros([2], dtype = np.int32)
        elem_v = xsp2xco1_v2type(vtype)
        xsp2co1_trace_98280(elem_v._data, tr, use_table)
        table[vtype] = tr[0]
        for j in range(2):
            w = rand_xsp2co1_elem()
            elem_v1 = w**-1 * elem_v  * w
            xsp2co1_trace_98280(elem_v1._data, tr[1:], use_table)
            assert tr[1] == tr[0]
    print(table)




#####################################################################
# Test conjugtion of extraspecial elements in group G_{x1}
#####################################################################

Co_1_ORDERS = set(list(range(1,17)) + 
      [18,20,21,22,23,24,26,28,30,33,35,36,39,40,42,60])
Gx0_ORDERS = set([x*y for x in Co_1_ORDERS for y in [1, 2, 4]])



