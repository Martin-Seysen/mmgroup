from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup import Xs12_Co1, PLoop, AutPL, Cocode
from mmgroup.generators import gen_leech2_type

#####################################################################
# Create test matrices
#####################################################################


def rand_xs1co1_elem():
    return Xs12_Co1(*[(x,) for x in "dxpylpylpylpy"])

def rand_xs1_vectors(length):
    return [randint(0, 0x1ffffff) for i in range(length)]

def create_conjugate_data():   
    for i in range(20):
        g = rand_xs1co1_elem()
        xs = rand_xs1_vectors(8)
        yield g, xs



#####################################################################
# Test conjugtion of extraspecial elements in group G_{x1}
#####################################################################

Co_1_ORDERS = set(list(range(1,17)) + 
      [18,20,21,22,23,24,26,28,30,33,35,36,39,40,42,60])
Gx0_ORDERS = set([x*y for x in Co_1_ORDERS for y in [1, 2, 4]])

    
@pytest.mark.qstate
def test_xs1_conjugate(verbose = 0):
    """Test the conjugation of Pauli matrix with unitary matrix"""
    for ntest, (g, xs) in enumerate(create_conjugate_data()):
        xs_g_all = [Xs12_Co1.from_xsp(x) for x in xs]
        xs_all = [x.as_xsp() for x in xs_g_all]
        ok = xs == xs_all
        o = g.order()
        assert o in Gx0_ORDERS
        if verbose or not ok:
            print("\nTest %d:\ng = %s" % (ntest, g))
            print("g has order %d" % o)
            print("v =", [hex(x) for x in xs])
            if not ok:
                print("Input and recomputed xs1 vectors v:")
                print("v =", [hex(x) for x in xs_all])
                err = "Error in recomputation of extraspecial vector"
                raise ValueError(err)

        conj = g.xsp_conjugate(xs)
        ref_conj = [(g**-1 * x * g).as_xsp()  for x in xs_g_all]
        ok = conj == ref_conj
        if verbose or not ok:
            if not verbose:
                print("Conjugation of v with g =\n", g)
            print("v =", [hex(x) for x in xs])
            print("Conjugated expected:\n", [hex(x) for x in ref_conj])
            if not ok:
                print("Conjugated obtained:\n", [hex(x) for x in conj])
                err = "Error in conjugation in group G_{x1}"
                raise ValueError(err)


#####################################################################
# Test computation of type Leech lattice vector
#####################################################################

ZERO = PLoop([])
OMEGA = ~ZERO
DODECAD = PLoop([0,4,8, 13,14,15, 17,18,19, 21,22,23])
OCTAD = PLoop([0,1,2,3,4,5,6,7])
HEXADECAD = ~OCTAD

def rand_n_elem():
    return Xs12_Co1(*[(x,) for x in "dxpy"])


def xs_vector(pl, cocode):
    return Xs12_Co1(("x", pl), ("d", Cocode(cocode))).as_xsp() 

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
    print ("Cocode:", cocode.syndrome_list(pos).bit_list)

def check_leech_type(x, t_expected):
    t = gen_leech2_type(x)
    ok = t == t_expected
    if not ok:
        print("Error: expected Leech type: %s, obtained: %s" % (
            hex(t_expected), hex(t)))
        display_leech_vector(x)
        err = "Error in computing Leech type"
        raise ValueError(err)
        

@pytest.mark.qstate
def test_xs1_type(verbose = 0):
    for ntest, (pl, cocode, vtype) in enumerate(type_data):
        x = xs_vector(pl, cocode)
        if verbose:
            print("\nTest %d", ntest+1)
            print("Expected type", hex(vtype))
            display_leech_vector(x)
        check_leech_type(x, vtype)
        for i in range(100):
            g =  rand_n_elem()
            x = g.xsp_conjugate(x)
            if verbose:
                display_leech_vector(x)
            check_leech_type(x, vtype)
