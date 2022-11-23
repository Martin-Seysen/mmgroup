from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral

import numpy as np
import pytest

from mmgroup import Xsp2_Co1, XLeech2
from mmgroup.clifford12 import xsp2co1_Co1_get_mapping
from mmgroup.clifford12 import xsp2co1_Co1_matrix_to_word
from mmgroup.clifford12 import xsp2co1_elem_from_mapping





#####################################################################
# Creating test vectors
#####################################################################


STD_BASIS = [1 << i  for i in range(24)]

PLUS = Xsp2_Co1()
MINUS = -PLUS
PM = [PLUS, MINUS]



def create_test_bases_get_mapping():
    def ar(m):
        return np.array(m, dtype = np.uint32)
    def make_from(c, t):
        m =  c.xsp_conjugate(STD_BASIS)
        m1 = t.xsp_conjugate(STD_BASIS)
        m2 = c.xsp_conjugate(m1)
        #assert m2 == t.c.xsp_conjugate(m)
        return ar(m1), ar(m2), ar(m), c


    testdata = [
        [[1 << i for i in range(24)]]*3 + [Xsp2_Co1()]
    ]
    for m1, m2, m, g in testdata:
        yield ar(m1), ar(m2), ar(m), g
    for i in range(100):
        c = Xsp2_Co1('r', 'G_x0') 
        t = Xsp2_Co1('r', 'G_x0')
        yield make_from(c, t)
 


def elem_from_mapping(m1, m2):
    a = np.zeros(10, dtype = np.uint32)
    res = xsp2co1_elem_from_mapping(m1, m2, a)
    if res < 0:
        err = "xsp2co1_elem_from_mapping returns %d"
        raise ValueError(err % res)
    g = Xsp2_Co1('a', a[:res & 0xff])
    order = (res >> 8) & 0xff 
    o = order //  (order & -order)
    assert g.order() == order
    nonzero = bool((res >> 16) & 1)
    chars = (g**o).chi_G_x0()
    chi = chars[2] * chars[3]
    assert chi >= 0
    assert nonzero == bool(chi)
    return g      




@pytest.mark.xsp2co1
def test_xsp2co1_Co1_get_mapping(verbose = 0):
    for n, (m1, m2, mref, g) in enumerate(create_test_bases_get_mapping()):
        m = np.zeros(24, dtype = np.uint32)
        res = xsp2co1_Co1_get_mapping(m1, m2, m)

        ok = res >= 0 and (m == mref).all()
        if verbose or not ok:
            print("Test", n)
            print('m1='   , m1)
            print('m2='   , m2)
            print('m_obt=', m)
            print('m_exp=', mref)
            if not ok:
                s = "Function xsp2co1_Co1_get_mapping returns %d"
                print(s % res)
                ERR = "Error in conjugating with element l"
                raise ValueError(ERR)
        o2, mod = divmod(g.order(), res)
        assert mod == 0, (g.order(), res)
        assert o2 > 0 and (o2 & (o2 - 1)) == 0, (g.order(), res)

        assert elem_from_mapping(m1, m2) * g**(-1) in PM
        


  


def create_test_elements(n = 100):
    def elem_to_matrix(elem):
        m = elem.xsp_conjugate([1 << i for i in range(24)])
        return np.array(m, dtype = np.uint32)
         
    testdata = [
    ]
    for elem in testdata:
        elem1 = Xsp2_Co1(elem)
        yield elem1, elem_to_matrix(elem1)
    for i in range(n):
        elem = Xsp2_Co1('r', 'G_x0') 
        yield elem, elem_to_matrix(elem)
 


@pytest.mark.xsp2co1
def test_xsp2co1_matrix_to_word(verbose = 0):
    for n, (elem, m) in enumerate(create_test_elements()):
        g = np.zeros(10, dtype = np.uint32)
        res = xsp2co1_Co1_matrix_to_word(m, g)
        elem_obt = Xsp2_Co1('a', g[:res]) if res >= 0 else None
        ok = res >= 0 and elem_obt * elem**(-1) in PM
        if verbose or not ok:
            print("Test", n)
            print('elem=    '   , elem)
            print('map=     '   , m)
            print('obtained='   , elem_obt )
            if elem_obt:
                print('prod=    '   , elem * elem_obt**(-1) )
            if not ok:
                if (res < 0): 
                    s = "Function xsp2co1_Co1_matrix_to_word returns %d"
                    print(s % res)
                ERR = "Error in conjugating with element l"
                try:
                    from mmgroup.clifford12 import xsp2co1_Co1_debug_pool_mapping
                    pool = np.zeros(72, dtype = np.uint64)
                    ln = xsp2co1_Co1_debug_pool_mapping(pool, len(pool))
                    if ln:
                        print("Error pool")
                        for i in range(len(pool)):
                            print("%2d: 0x%08x" % (i,pool[i]))
                except:
                    pass
                raise ValueError(ERR)




    

    

