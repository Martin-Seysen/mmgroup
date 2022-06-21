from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from multiprocessing import Pool
import time

import numpy as np
import pytest

from mmgroup import Xsp2_Co1, PLoop, AutPL, Cocode, MM0, MM, XLeech2
from mmgroup.generators import gen_leech2_op_word
from mmgroup.clifford12 import xsp2co1_set_elem_word_scan
from mmgroup.clifford12 import xsp2co1_unit_elem
from mmgroup.structures.construct_mm import iter_mm


#####################################################################
# Create test matrices
#####################################################################



def construct_xsp2co1_testcases():
    def a(string_data):
        data = [(x, 'n') for x in string_data]
        return np.array(list(iter_mm(None, data)), dtype = np.uint32)

    for tag0 in "pdxyl":
        yield a(tag0), 1
        for tag1 in "pdxyl":
            yield a(tag0 + tag1 + 't'), 2

    for i in range(5):
        data = Xsp2_Co1('r', 'G_x0').mmdata
        yield data, len(data)

def dirty_elem():
    return np.array([randint(0, 0xffffffffffffffff) for i in range(26)],
        dtype = np.uint64)

@pytest.mark.xsp2co1
def test_xsp2co1_set_elem_word_scan():
    for a, len_a in construct_xsp2co1_testcases():
        elem0 = dirty_elem()
        elem1 = dirty_elem()
        s0 =  xsp2co1_set_elem_word_scan(elem0, a, len(a), 0)
        assert s0 == len_a
        xsp2co1_unit_elem(elem1)
        s1 =  xsp2co1_set_elem_word_scan(elem1, a, len(a), 1)
        assert s1 == len_a
        assert s0 == s1
        xsp = Xsp2_Co1()
        np.copyto(xsp._data, elem0)
        assert Xsp2_Co1('a', a[:len_a]) == xsp






