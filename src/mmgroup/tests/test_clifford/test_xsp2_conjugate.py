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
from mmgroup.generators import gen_leech2_op_word_leech2
from mmgroup.generators import gen_leech2_type
from mmgroup.clifford12 import xsp2co1_isotropic_type4
from mmgroup.clifford12 import bitmatrix64_vmul


###################################################################
###################################################################
# Test C function gen_leech2_op_word and derived python functions
###################################################################
###################################################################


#####################################################################
# Create test matrices
#####################################################################


def rand_xsp2co1_elem():
    return MM0([(x,'r') for x in "dxpylpylpylpy"])

def rand_xsp2_vectors(length):
    return [randint(0, 0x1ffffff) for i in range(length)]

def create_conjugate_data():  
    # Testing with following test data verfies that generators
    # :math:`\xi^e, e = 1,2` operate on the group :math:`Q_{x0}`
    # correctly. This verifies the correct implementation of 
    # these generator up to sign. The sign can be verified by
    # checking that the generator has order 3
    for exp in (1, 2): 
        g_mm = MM0('l', exp)
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
    w_g = g_mm.mmdata
    return gen_leech2_op_word(x, w_g, len(w_g))


    
@pytest.mark.xsp2co1
def test_xsp2_conjugate(verbose = 0):
    """Test the conjugation of Pauli matrix with unitary matrix"""
    l0, l1, l2 = Xsp2_Co1(), Xsp2_Co1('l', 1),  Xsp2_Co1('l', 2)
    assert l1**2 == l2
    assert l1**3 == l0    
    for ntest, (g_mm, xs) in enumerate(create_conjugate_data()):
        g = Xsp2_Co1(g_mm)
        xs_g_all = [Xsp2_Co1.group.from_xsp(x) for x in xs]
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


###################################################################
###################################################################
# Test C function gen_leech2_op_word_Omega
###################################################################
###################################################################


def gen_leech2_op_word_Omega_testdata():  
    Omega = 0x800000
    testdata = [[('l', i)] for i in range(1,2)]
    tags = "xyzdp"
    testdata += [[(tag, 'n')] for tag in tags]
    for data in testdata:
        yield Omega, MM0(data)
    testdata2 = [[('l', 'n'),(tag, 'n'),('l', 'n')] for tag in tags]
    for data in testdata2:
        g = MM0(data)
        yield Omega, g
        a = g.mmdata
        a[1] ^= 0x80000000
        g = MM0('a', a)
        yield Omega, g
    # test with some more random data
    for i in range(10):
        g = rand_xsp2co1_elem()
        yield Omega, g
        yield randint(0, 0xffffff), g


@pytest.mark.xsp2co1
def test_gen_leech2_op_word_leech2(verbose = 0):
    for ntest, (l, g) in enumerate(gen_leech2_op_word_Omega_testdata()):
        for inv in [0, 1]:
            g1 = g ** -1 if inv else g
            gdata, len_gdata = g.mmdata, len(g.mmdata)
            ref = conj_x_by_word(l, g1 ) & 0xffffff
            res = gen_leech2_op_word_leech2(l, gdata, len(gdata), inv)
            ok = ref == res
            if verbose or not ok: 
                inv = " (inverse)" if inv else ""
                print("\nTest %d%s:\n g = %s" % (ntest, inv, g))
                print(" Leech2 vector", hex(l))
                if len(g.mmdata) <= 3:
                    print(" Data of g:", [hex(x) for x in g1.mmdata])
                if ok:
                    print(" Result:", hex(res))
                else:
                    print(" Obtained:", hex(res))
                    print(" Expected:", hex(ref))
                    err = "Error in C function gen_leech2_op_word_Omega"
                    raise ValueError(err)

    



###################################################################
###################################################################
# xsp2co1_isotropic_type4
###################################################################
###################################################################

def xsp2co1_isotropic_type4_testdata():
    """Testdata for function xsp2co1_isotropic_type4

    Yield pairs (L, n), where L is a list list of type-4 vectors,
    and `n` is the dimension of the result of function
    ``xsp2co1_isotropic_type4`` applied to that list.
    """
    Omega = XLeech2(0x800000, 0)
    omega =  XLeech2(0, [0,1,2,3])
    octad = XLeech2(PLoop(range(8))) * Omega
    print(Omega, omega)
    for i in range(30):
        g = Xsp2_Co1('r', 'G_x0')
        yield [Omega*g, omega*g], 8
    for i in range(30):
        g = Xsp2_Co1('r', 'G_x0')
        yield [Omega*g, octad*g], 6




def check_isotropic_type4(v, a):
    for w in v:
        v_data = v[0].ord
        for i in range(30):
            r = randint(0, (1 << len(v)) - 1)
            x = bitmatrix64_vmul(r, a, len(a))
            x_data = x
            w_data = w.ord
            t_x = gen_leech2_type(w_data)
            assert t_x & 1 == 0
            t_wx = gen_leech2_type(w_data ^ x_data)
            mask = 3 if len(v) >= 8 else 1
            assert (t_x ^ t_wx) & mask == 0

@pytest.mark.xsp2co1
def test_xsp2co1_isotropic_type4(verbose = 0):
    if verbose:
        print("Testing function xsp2co1_isotropic_type4")
    a = np.zeros(24, dtype = np.uint64)
    for v, n in xsp2co1_isotropic_type4_testdata():
        dim = -1
        if verbose:
            print("v=", v, "n =", n)
        assert  gen_leech2_type(v[0].ord) == 4
        dim = xsp2co1_isotropic_type4(v[0].ord, a, -1)
        assert dim == 12, dim
        check_isotropic_type4(v[:1], a[:dim])
        for w in v[1:]:
             assert  gen_leech2_type(w.ord) == 4
             dim = xsp2co1_isotropic_type4(v[1].ord, a, dim)
             pass
        assert dim == n
        check_isotropic_type4(v, a[:dim])

