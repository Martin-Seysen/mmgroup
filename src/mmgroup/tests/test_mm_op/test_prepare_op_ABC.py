from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from random import randint
from mmgroup import MM0, Xsp2_Co1, AutPL, MM

import pytest


from mmgroup.mm_op import mm_group_prepare_op_ABC
from mmgroup.generators import mm_group_n_clear
from mmgroup.generators import mm_group_n_mul_element
from mmgroup.generators import mm_group_n_mul_word_scan
from mmgroup.generators import mm_group_n_to_word
from mmgroup.generators import mm_group_n_right_coset_N_x0

from mmgroup.clifford12 import xsp2co1_set_elem_word_scan
from mmgroup.clifford12 import xsp2co1_elem_subtype
from mmgroup.clifford12 import xsp2co1_elem_to_N0
from mmgroup.clifford12 import xsp2co1_mul_elem_word
from mmgroup.clifford12 import xsp2co1_elem_to_word


def display_elem(elem):
    g = Xsp2_Co1()
    np.copyto(g._data, elem)
    print(g)
    return g
    


def py_prepare_op_ABC(g):
    a = np.zeros(11, dtype = np.uint32)

    g = g.mmdata
    assert g.dtype == np.uint32
    has_t = has_l = reduce = False
    for x in g:
        reduce |= has_t
        tag, data = (x >> 28) & 7, x & 0xfffffff
        if tag == 5 and data: 
            has_t = True
        if tag == 6 and data: 
            has_l = True
        if tag == 7 and has_data:
             raise ValueError("Illegal tag")

    if not has_l:
        mm_group_n_clear(a)
        len_a = mm_group_n_mul_word_scan(a, g, len(g))
        assert len_a == len(g) 
        len_a = mm_group_n_to_word(a, a)
        return 1, a[:len_a]

    if reduce or len(g) > 11:
        gn = np.zeros(5, dtype = np.uint32)
        elem = np.zeros(26, dtype = np.uint64)
        pos = xsp2co1_set_elem_word_scan(elem, g, len(g), 0)
        assert 0 <= pos <= len(g)
        #display_elem(elem)
        mm_group_n_clear(gn)
        scan = mm_group_n_mul_word_scan(gn, g[pos:], len(g) - pos)
        assert scan >= 0
        assert pos + scan == len(g)
        
        if xsp2co1_elem_subtype(elem) == 0x48:
            res = xsp2co1_elem_to_N0(elem, a)
            assert res >= 0
            mm_group_n_mul_element(a, gn, a)
            len_a = mm_group_n_to_word(a, a)
            return 1, a[:len_a]
        else:
            e = mm_group_n_right_coset_N_x0(gn)
            len_a = mm_group_n_to_word(gn, gn)
            res = xsp2co1_mul_elem_word(elem, gn, len_a)
            assert res >= 0
            len_a = xsp2co1_elem_to_word(elem, a)
            assert 0 <= len_a <= 10
            if e:
                a[len_a] = 0x50000000 + e
                len_a += 1
            return 0, a[:len_a]

    return 0, np.copy(g)
            



def prep_op_ABC(g, must_be_in_N0 = False, verbose = 1):
    a = np.zeros(11, dtype = np.uint32)
    py_prepare_op_ABC(g)
    res = mm_group_prepare_op_ABC(g.mmdata, len(g.mmdata), a)
    if res < 0:
        print("Test function mm_group_prepare_op_ABC with g =")
        print(g)
        s = "Function mm_group_prepare_op_ABC returns with error %d" 
        raise ValueError(s % res)
    b0, len_a =  (res >> 8) & 1, res & 0xff
    a = a[:len_a]
    g1 = type(g)('a', a)
    return b0, a


def prep_op_ABC_testdata():
    for i in range(20):
        yield MM0('r', 'G_x0') * MM0('r', 'N_0'), False
    for i in range(20):
        yield MM0('r', 'N_0'), True

    pi0 = zip(list(range(10)), [x ^ 1 for x in range(10)])
    #print(pi0)
    g =  MM0([('l', 1), ('p', pi0), ('d', 0x800), ('l', 2), ('d', 0x800)])**5
    assert len(g.mmdata) > 11
    _ = Xsp2_Co1((MM0(g).reduce()))  # This should exist!!
    yield g, True
    


@pytest.mark.mm_op
def test_prepare_op_ABC(verbose = 0):
    print("Testing function mm_group_prepare_op_ABC")
    for i, (g, must_be_in_N0) in enumerate(prep_op_ABC_testdata()):
        b0, g1 = py_prepare_op_ABC(g)
        try:
            ok = g == g.group('a', g1)
        except:
            print("Could not compare elements in Monster group")
            ok = False
        if must_be_in_N0:
            ok &= b0 == 1
        b02, g2 = prep_op_ABC(g, must_be_in_N0, verbose = verbose)
        #print(type(g2), type(g2), g1, g2)
        ok &= b02 == b0 and len(g1) == len(g2) and (g2 == g1).all()
        if  verbose or not ok:
            print("Input g:\n" + str(g))
            print("Py Output g (b0 = %d):\n" % b0 + str(g1))
            if must_be_in_N0 and not b0:
                print("b0 should be set but isn't")
            if (g2 != g1).any() or b02 != b0 :
                print("C Output g (b0 = %d):\n" % b02 + str(g2))
            if not ok:
                raise ValueError("Function py_prepare_op_ABC failed")
        b02, g2 = prep_op_ABC(g, must_be_in_N0, verbose = verbose)
       
    print("Test passed")





