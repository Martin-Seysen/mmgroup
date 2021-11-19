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
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_type2
from mmgroup.generators import gen_leech2_count_type2
from mmgroup.generators import gen_leech2_op_word
from mmgroup.clifford12 import xsp2co1_leech2_count_type2



#####################################################################
# Test computation of type Leech lattice vector
#####################################################################

ZERO = PLoop([])
OMEGA = ~ZERO
DODECAD = PLoop([0,4,8, 13,14,15, 17,18,19, 21,22,23])
OCTAD = PLoop([0,1,2,3,4,5,6,7])
HEXADECAD = ~OCTAD

def rand_n_elem():
    return Xsp2_Co1([(x,'r') for x in "dxpy"])


def xs_vector(pl, cocode):
    return Xsp2_Co1([("x", pl), ("d", Cocode(cocode))]).as_xsp() 

type_data = [
    (ZERO, [], 0),
    (OMEGA, [], 0x48),
    (ZERO, [3,4], 0x20),
    (ZERO, [1], 0x21),
    (OCTAD, [], 0x22),
    (OMEGA, [1], 0x31),
    (OCTAD, [7,8], 0x34),
    (OMEGA, [1,2,3], 0x33),
    (DODECAD, [], 0x36),
    (ZERO, [0,1,2,3], 0x40),
    (OCTAD, [0,1], 0x42),
    (ZERO, [1,2,3], 0x43),
    (OCTAD, [8,9], 0x44),
    (DODECAD, [0,1], 0x46),
]
   


def display_leech_vector(x):
    gcode = PLoop(x >> 12)
    bl = gcode.bit_list
    print ("GCode:\n", bl)
    if len(bl) == 16:
        l = [x for x in range(24) if not x in bl]
        pos = bl[0]
    elif len(gcode):
        pos = bl[0]
    else:
        pos = 0
    cocode = Cocode(x) + gcode.theta()
    print ("Cocode:", cocode.syndrome_list(pos))


def alternative_is_type2(v):
    a = np.array([0, v], dtype = np.uint64)
    return xsp2co1_leech2_count_type2(a, 2)

def check_leech_type(x, t_expected):
    t = gen_leech2_type(x)
    ok = t == t_expected
    if not ok:
        print("Error: expected Leech type: %s, obtained: %s" % (
            hex(t_expected), hex(t)))
        display_leech_vector(x)
        err = "Error in computing Leech type"
        raise ValueError(err)
    is_type2 = (t_expected >> 4) == 2
    found_type2 = gen_leech2_type2(x) > 0
    ok = is_type2 == found_type2
    if not ok:
        print("Error:  x = %s, Leech type: %s" % (
            hex(x),  hex(t_expected)), is_type2, hex(found_type2))
        err = "Function gen_leech2_type2 failed"
        raise ValueError(err)
    alt_found_type2 = alternative_is_type2(x)
    ok = is_type2 == alt_found_type2
    if not ok:
        print("Error:  x = %s, Leech type: %s" % (
            hex(x),  hex(t_expected)), is_type2, alt_found_type2)
        err = "Function xsp2co1_leech2_count_type2 failed"
        #raise ValueError(err)
        

@pytest.mark.xsp2co1
def test_xsp2_type(verbose = 0):
    for ntest, (pl, cocode, vtype) in enumerate(type_data):
        x = xs_vector(pl, cocode)
        if verbose:
            print("\nTest %d" %(ntest+1))
            print("Expected type", hex(vtype))
            display_leech_vector(x)
        check_leech_type(x, vtype)
        for i in range(200):
            g =  rand_n_elem()
            x = g.xsp_conjugate(x)
            if verbose:
                display_leech_vector(x)
            check_leech_type(x, vtype)




#####################################################################
# Test function gen_leech2_type2 via selftest in C file
#####################################################################


@pytest.mark.xsp2co1
def test_leech2_self(fast = 1, verbose = 0):
    f = xsp2co1_leech2_count_type2 if fast else gen_leech2_count_type2
    a_type = np.uint64 if fast else np.uint32
    if verbose:
        print("Testing gen_leech2_type2() ... ")
    a = np.array([0]+[1 << i for i in range(24)], dtype = a_type)
    t_start = time.time()
    result = f(a, 25)
    t = time.time() - t_start
    assert result == 98280, result
    if verbose:
        print("passed, %.3f ms" % (1000*t))

#####################################################################
# Test method subtype of class Xsp2Co1
#####################################################################


def subtype_testdata():
    yield XLeech2(~PLoop())
    yield XLeech2(Cocode([0,1,2,3]))
    yield XLeech2(~PLoop(list(range(8))))
    yield XLeech2(~PLoop(list(range(8))), Cocode([8,9]))
    for i in range(50):
        yield XLeech2('r', 4)


@pytest.mark.xsp2co1
def test_subtype():
    OMEGA = XLeech2(~PLoop())
    types = set()
    for v in subtype_testdata():
         g = MM('c', v)
         v2 =  Xsp2_Co1(g).subtype
         v2ref = (OMEGA * g).subtype
         assert v2 == v2ref, (v2, v2ref)
         types.add(v2)
    assert len(types) == 6



