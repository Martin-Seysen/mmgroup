from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from multiprocessing import Pool

import numpy as np
import pytest

from mmgroup import Xsp2_Co1, PLoop, AutPL, Cocode, MM
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_is_type2
from mmgroup.generators import gen_leech2_is_type2_selftest
from mmgroup.generators import gen_leech2_op_word



#####################################################################
# Create test matrices
#####################################################################


def rand_xsp2co1_elem():
    return MM(*[(x,) for x in "dxpylpylpylpy"])

def rand_xsp2_vectors(length):
    return [randint(0, 0x1ffffff) for i in range(length)]

def create_conjugate_data():  
    # Testing with following test data verfies that generators
    # :math:`\xi^e, e = 1,2` operate on the group :math:`Q_{x0}`
    # correctly. This verifies the correct implementation of 
    # these generator up to sign. The sign can be verified by
    # checking that the generator has order 3
    for exp in (1, 2): 
        g_mm = MM(('l', exp))
        xs = [1 << i for i in range(24)]
        yield g_mm, xs         
    # test with some more random data
    for i in range(20):
        g = rand_xsp2co1_elem()
        xs = rand_xsp2_vectors(8)
        yield g, xs



#####################################################################
# Test conjugtion of extraspecial elements in group G_{x1}
#####################################################################

Co_1_ORDERS = set(list(range(1,17)) + 
      [18,20,21,22,23,24,26,28,30,33,35,36,39,40,42,60])
Gx0_ORDERS = set([x*y for x in Co_1_ORDERS for y in [1, 2, 4]])


def conj_x_by_word(x, g_mm):
    """Conugate x ith a word of elments in the Clifford group

    Here :math:`x` is an element of the normal subgroup :math:`Q_x0`
    of :math:`G_x0`, and :math:`g_{mm}` is an element of the 
    subgroup :math:`G_x0` of the monster.

    :math:`x` is given as an integer in **Leech lattice
    represention** and :math:`g_{mm}` is given as an element
    of the group  object ``mmgroup.MM``, which is the standard 
    instance of the monster group.

    The function returns :math:`g_{mm}^{-1} x g_{mm}` as an
    integer in **Leech lattice represention**.
    """
    w_g = g_mm.data
    return gen_leech2_op_word(x, w_g, len(w_g))

    
@pytest.mark.qstate
def test_xsp2_conjugate(verbose = 0):
    """Test the conjugation of Pauli matrix with unitary matrix"""
    l0, l1, l2 = Xsp2_Co1(), Xsp2_Co1(('l', 1)),  Xsp2_Co1(('l', 2))
    assert l1**2 == l2
    assert l1**3 == l0    
    for ntest, (g_mm, xs) in enumerate(create_conjugate_data()):
        g = Xsp2_Co1(g_mm)
        xs_g_all = [Xsp2_Co1.from_xsp(x) for x in xs]
        xs_all = [x.as_xsp() for x in xs_g_all]
        ok = xs == xs_all
        o = g.order()
        assert o in Gx0_ORDERS
        if verbose or not ok:
            print("\nTest %d:\ng = %s" % (ntest, g))
            print("g has order %d" % o)
            print("v =", [hex(x) for x in xs])
            if not ok:
                print("Input and recomputed xsp2 vectors v:")
                print("v =", [hex(x) for x in xs_all])
                err = "Error in recomputation of extraspecial vector"
                raise ValueError(err)

        conj = g.xsp_conjugate(xs)
        ref_conj = [(g**-1 * x * g).as_xsp()  for x in xs_g_all]
        ref_conj1 = [conj_x_by_word(x, g_mm) for x in xs]
        ref_ok = ref_conj == ref_conj1
        if not ref_ok:
            print("Conjugation of v with g =\n", g)
            print("v =", [hex(x) for x in xs])
            print("Reference for v * g via direct calculation")
            print([hex(x) for x in ref_conj]) 
            print("Reference using C function gen_leech2_op_word")
            print([hex(x) for x in ref_conj1]) 
            err = "Error in reference for conjugation in group G_{x1}"
            raise ValueError(err)

        ok = ref_ok and conj == ref_conj1

        if verbose or not ok:
            if not verbose:
                print("Conjugation of v with g =\n", g)
            print("v =", [hex(x) for x in xs])
            print("Conjugated expected:\n", [hex(x) for x in ref_conj])
            if not ok:
                print("Conjugated obtained:\n", [hex(x) for x in conj])
                err = "Error in conjugation in group G_{x1}"
                raise ValueError(err)

        conj_nosign = g.xsp_conjugate(xs, sign = False)
        mask = (1 << 24) - 1
        assert conj_nosign == [x & mask for x in conj]


#####################################################################
# Test computation of type Leech lattice vector
#####################################################################

ZERO = PLoop([])
OMEGA = ~ZERO
DODECAD = PLoop([0,4,8, 13,14,15, 17,18,19, 21,22,23])
OCTAD = PLoop([0,1,2,3,4,5,6,7])
HEXADECAD = ~OCTAD

def rand_n_elem():
    return Xsp2_Co1(*[(x,) for x in "dxpy"])


def xs_vector(pl, cocode):
    return Xsp2_Co1(("x", pl), ("d", Cocode(cocode))).as_xsp() 

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
    is_type2 = (t_expected >> 4) == 2
    found_type2 = gen_leech2_is_type2(x) 
    ok = is_type2 == found_type2
    if not ok:
        print("Error:  v = %s, Leech type: %s" % (
            hex(x),  hex(t_expected)), is_type2, hex(found_type2))
        err = "Function gen_leech2_is_type2 failed"
        raise ValueError(err)
        

@pytest.mark.qstate
def test_xsp2_type(verbose = 0):
    for ntest, (pl, cocode, vtype) in enumerate(type_data):
        x = xs_vector(pl, cocode)
        if verbose:
            print("\nTest %d", ntest+1)
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
# Test function gen_leech2_is_type2 via selftest in C file
#####################################################################

def one_selftest_leech_type2(data):
    return  gen_leech2_is_type2_selftest(*data)


def gen_selftest_inputs(n):
    assert 0x1000000 % n == 0
    q = 0x1000000 // n
    for  i in range(n):
        yield i*q, q


@pytest.mark.qstate
def test_leech2_self(verbose = 0):
    NPROCESSES = 4
    if verbose:
        print("Testing gen_leech2_is_type2() ... ", end = "")
    with Pool(processes = NPROCESSES) as pool:
        results = pool.map(one_selftest_leech_type2, 
                   gen_selftest_inputs(NPROCESSES))
    pool.join()
    result = sum(results)
    assert result == 98280, result
    if verbose:
        print("passed")
