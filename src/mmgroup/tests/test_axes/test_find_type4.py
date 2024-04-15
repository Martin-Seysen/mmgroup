
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import sample, randint, shuffle
from collections import defaultdict

import numpy as np
import pytest

from mmgroup import MM, MMV, Xsp2_Co1, XLeech2
from mmgroup.generators import gen_leech2_subtype, gen_leech2_type


def import_all():
    global mm_reduce_find_type4
    from mmgroup.mm_reduce import mm_reduce_find_type4


VALUE_T4 = defaultdict(lambda : 5)
VALUE_T4.update({0x48:0, 0x40:1, 0x42:2, 0x44:2, 0x46:3, 0x43:4})
#print(VALUE_T4)


def eval_type4(v, v2 = 0):
    value = VALUE_T4[gen_leech2_subtype(v)]
    if v2 & 0xffffff and gen_leech2_type(v ^ v2) != 2:
        value = 5
    return min(0x5000000, (v & 0xffffff) + (value << 24))


#####################################################################################
# Test data for function mm_reduce_find_type4
#####################################################################################

def rand_type_40_vector():
    return randint(0, 0xffffff) & 0x8007ff


def rand_type2_type4_list(length, random_length = 0):
    v2 = XLeech2('r', 2)
    v4_list = []
    while(len(v4_list) < length):
       v4 = XLeech2('r', 2) * v2
       if v4.type == 4:
           v4_list.append(v4.ord)
    v4_list += [randint(0, 0xffffff) for x in range(random_length)]
    shuffle(v4_list)
    return v2.ord, np.array(v4_list, dtype = np.uint32)
   


def find_type4_testcases(ncases = 1):
    data = [
        ([0x800000, 0x800700, 0], 0,  0x800000), 
    ]
    for a, v2, v_best in data:
        yield np.array(a, dtype = np.uint32), v2, v_best
    for i in range(500):
        a = [rand_type_40_vector()] + [randint(0x1000, 0x7fffff) 
            for i in range(4)]
        best = a[0] if gen_leech2_type(a[0]) == 4 else None
        shuffle(a)
        yield np.array(a, dtype = np.uint32), 0, best
    lengths = list(range(20, 100, 3)) + list(range(2000, 2100,7))
    for i in lengths:
        a = np.random.randint(0, 0xffffff, i, dtype = np.uint32)
        yield a, 0, None
    for i in range(1000):
        v2, a = rand_type2_type4_list(3, 5)
        yield a, v2, None
             


#####################################################################################
# Test function mm_reduce_find_type4
#####################################################################################


@pytest.mark.mmgroup
def test_find_type4(verbose = 0):
    import_all()
    print("Test function mm_reduce_find_type4")
    for a, v2, v_best  in find_type4_testcases():
        if verbose:
            print("v2 =", hex(v2), ", v_best =", v_best)
            if len(a) <= 10: print("a=", a)
        v0 =  mm_reduce_find_type4(a, len(a), v2)
        val_0 = eval_type4(v0, v2)
        val_a = [eval_type4(x, v2) for x in a]
        if verbose:
            if len(a) <= 10: 
                print("a=", a)
                print("subtypes =", [hex(gen_leech2_subtype(x)) for x in a])
                if v2 & 0xffffff:
                    print("subtypes_2 =", [
                        hex(gen_leech2_subtype(x ^ v2)) for x in a])
            print("best =", hex(v0), ", val_0=", hex(val_0))
            #print("val_a =", [hex(x) for x in val_a])
        assert val_0 <= min(val_a), (hex(val_0), hex(min(val_a)))
        if v_best is None:
            if val_0 >> 24 == 5:
                assert v0 == 0, (hex(v0), gen_leech2_subtype(v0))
            else:
                assert v0 in a, hex(v0)
        else:
            assert v_best == v0, (hex(v_best), hex(v0))





        

