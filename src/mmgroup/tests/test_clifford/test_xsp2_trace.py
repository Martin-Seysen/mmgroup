from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from multiprocessing import Pool
from collections import defaultdict
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
from mmgroup.clifford12 import xsp2co1_trace_4096
from mmgroup.clifford12 import xsp2co1_traces_small
from mmgroup.clifford12 import xsp2co1_traces_all

from mmgroup.tests.test_mm.check_monster_orders import ClassOrders
from mmgroup.tests.test_mm.check_monster_orders import CharacterValues

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


@pytest.mark.xsp2co1
def test_xsp2_count_table():
    use_table = 1
    table = [0] * 5
    for vtype in [0,2,3,4]:
        tr = np.zeros([2], dtype = np.int32)
        elem_v = xsp2xco1_v2type(vtype)
        xsp2co1_trace_98280(elem_v._data, tr, use_table)
        table[vtype] = tr[0]
        for j in range(2):
            w = rand_xsp2co1_elem()
            elem_v1 = w**-1 * elem_v  * w
            xsp2co1_trace_98280(elem_v1._data, tr[1:], 0)
            assert tr[1] == tr[0]
    print(table)




#####################################################################
# Test conjugtion of extraspecial elements in group G_{x1}
#####################################################################



assert len(ClassOrders) == len(CharacterValues)

CharacterDict = defaultdict(set)
for order, char_value in zip(ClassOrders, CharacterValues):
    CharacterDict[order].add(char_value)

#print(CharacterDict)
 

Xsp2_Co1_Element = type(Xsp2_Co1())

def character(g, verbose = 0):
    assert isinstance(g, Xsp2_Co1_Element)
    a = np.zeros(4, dtype = np.int32)
    res = xsp2co1_traces_all(g._data, a)
    if verbose: print("Subcharacters", [int(x) for x in a])
    assert res >= 0, res
    chi24, chisq24, chi4096, chi98260 = map(int, a[:4])
    chi299 = (chi24**2 + chisq24) // 2 - 1
    assert chi24 >= 0
    if chi24 == 0: assert chi4096 >= 0
    return chi299 + chi98260 + chi24 * chi4096


def character_testcases():
    data = [
        [],
        [("x", 0x1f4c), ("d", 0x375)],
    ]
    for d in data:
        yield Xsp2_Co1(*d)
    for i in range(100):
        yield  rand_xsp2co1_elem()

@pytest.mark.xsp2co1
def test_xsp2_characters(verbose = 0):
    if verbose:
        print("Test character calculation in G_x0")
    for n, g in enumerate(character_testcases()):
        if verbose: print("Test %d, g=" % (n+1), g)
        o = g.order()
        good_orders = CharacterDict[o] 
        if verbose: print("Order =", o, ", chars =", good_orders)
        chi = character(g, verbose)
        if verbose:  print("Character =", chi)
        assert chi in good_orders, (o, chi, good_orders)



