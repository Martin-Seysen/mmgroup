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

from mmgroup import Xsp2_Co1, PLoop, AutPL, Cocode, MM0
from mmgroup.generators import gen_leech2_type
from mmgroup.clifford12 import xsp2co1_traces_all
from mmgroup.clifford12 import xsp2co1_traces_fast

from mmgroup.tests.test_orders.check_monster_orders import ClassOrders
from mmgroup.tests.test_orders.check_monster_orders import CharacterValues

from mmgroup.tests.test_involutions.test_involution_invariants import INVOLUTION_SAMPLES


#####################################################################
# Create test matrices
#####################################################################



def rand_xsp2co1_elem():
    return Xsp2_Co1('r', 'G_x0')

def xsp2xco1_xsp2(v):
    pl = PLoop(v >> 12)
    coc = Cocode(v)
    return Xsp2_Co1([("x", v >> 12), pl.theta(), coc])


v2_types = {0:0}


def xsp2xco1_v2type(vtype):
    global v2_types
    try:
        return xsp2xco1_xsp2(v2_types[vtype])
    except KeyError:
        assert vtype in [2,3,4] 
        for i in range(10000):
            v = randint(0, 0xffffff)
            if gen_leech2_type(v) == vtype:
                v2_types[vtype] = xsp2xco1_xsp2(v) 
                return v2_types[vtype]
        raise ValueError("No Leech lattice vector of type", vtype)



def xsp2co1_traces(elem):
    tr = np.zeros([4], dtype = np.int32)
    tr_fast = np.zeros([4], dtype = np.int32)
    res = xsp2co1_traces_all(elem._data, tr)
    assert res >= 0
    res_fast = xsp2co1_traces_fast(elem._data, tr_fast)
    assert res_fast >= 0
    assert (tr == tr_fast).all()
    return [int(x) for x in tr]



def xsp2co1_trace_98280(elem):
    return xsp2co1_traces(elem)[3]
     

@pytest.mark.involution
def test_xsp2_count_table():
    table = [0] * 5
    for vtype in [0,2,3,4]:
        elem_v = xsp2xco1_v2type(vtype)
        tr_ref = xsp2co1_trace_98280(elem_v)
        table[vtype] = tr_ref
        for j in range(2):
            w = rand_xsp2co1_elem()
            elem_v1 = w**-1 * elem_v  * w
            tr = xsp2co1_trace_98280(elem_v1)
            assert tr == tr_ref
    print(table)





#####################################################################
# Test computation of characters in group G_{x0}
#####################################################################



assert len(ClassOrders) == len(CharacterValues)

CharacterDict = defaultdict(set)
for order, char_value in zip(ClassOrders, CharacterValues):
    CharacterDict[order].add(char_value)

#print(CharacterDict)
 

Xsp2_Co1_Element = type(Xsp2_Co1())

def character(g, verbose = 0):
    assert isinstance(g, Xsp2_Co1_Element)
    a = xsp2co1_traces(g)
    if verbose: print("Subcharacters", a)
    chi24, chisq24, chi4096, chi98260 = map(int, a[:4])
    chi299 = (chi24**2 + chisq24) // 2 - 1
    assert chi24 >= 0
    if chi24 == 0: assert chi4096 >= 0
    chi_M = chi299 + chi98260 + chi24 * chi4096
    #print("MMMM", MM0(g))
    chi_g = MM0(g).chi_G_x0()
    assert chi_g == (chi_M, chi299, chi24, chi4096)
    return chi_M


def character_testcases():
    data = [
        [],
        [("x", 0x1f4c), ("d", 0x375)],
    ]
    for d in data:
        yield Xsp2_Co1(d)
    for i in range(200):
        yield  rand_xsp2co1_elem()

@pytest.mark.involution
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




@pytest.mark.involution
def test_involution_samples(verbose = 0):
    if verbose:
        print("Test character calculation on involutions in G_x0/Q_x0")
    for (_1, chi, _2), g_str in  INVOLUTION_SAMPLES:
        chi_576 = 2 * (chi[1]  + 1) - chi[2]**2
        chi_ref = [chi[2], chi_576, chi[3], chi[4]]
        g = Xsp2_Co1(g_str)
        for j in range(4):
            w = rand_xsp2co1_elem()
            assert xsp2co1_traces(g**w) == chi_ref


