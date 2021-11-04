from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from random import randint

import pytest


from mmgroup.mm import mm_sub_get_table_xi
from mmgroup.tests.groups.mgroup_n import MGroupNWord
from mmgroup.tests.spaces.sparse_mm_space import SparseMmV


space = SparseMmV(3)
group = MGroupNWord

 

def unit_vector_to_tuple(v):
    tuples = v.as_tuples()
    assert len(tuples) == 1
    t = tuples[0]
    assert t[0] == 1
    return t[1:]


def tuple_to_xi_tuple(t):
    tag, i, j = t
    if tag == 'B':
        return "BC", 32 * i + j
    if tag == 'C':
        return "BC", 32 * (i + 24) + j
    if tag == 'T':
        if 0 <= i < 15:
            return "BC", 32 * 48 + i *64 + j
        if 15 <= i < 375:
            return "T0", 64 * (i - 15) + j
        return "T1", 64 * (i - 375) + j
    if tag == 'X':
        if i < 1024:
            return "X0", 32 * i + j
        return "X1", 32 * (i - 1024) + j


def xi_tuple_to_tuple(t):
    tag, i = t
    if tag == "BC":
        if i >= 48 * 32:
            i1 = i - 48 * 32
            return "T", i1 >> 6, i1 & 63
        if i > 24 * 32:
            i1 = i - 24 * 32
            return "C", i1 >> 5, i1 & 31
        return "B", i >> 5, i & 31
    if tag == "T0":
        return "T", (i >> 6) + 15, i & 63
    if tag == "T1":
        return "T", (i >> 6) + 375, i & 63
    if tag == "X0":
        return "X", (i >> 5), i & 31
    if tag == "X1":
        return "X", (i >> 5) + 1024, i & 31


def op_xi_tuple(v, e, verbose =  0):
    from mmgroup.dev.mm_basics.mm_tables_xi import MM_TablesXi
    tables_xi = MM_TablesXi()
    e1 = e - 1
    table_index = None
    for i, source_tag in tables_xi.SOURCE_TAGS.items():
        if source_tag[e1] == v[0]:
            table_index = i
    dest_tag = tables_xi.DEST_TAGS[table_index][e1]
    dest_shape = tables_xi.DEST_SHAPES[table_index]
    source_shape = tables_xi.SOURCE_SHAPES[table_index]
    block, entry = divmod(v[1], source_shape[1] * 32)
    section_len = dest_shape[1] * dest_shape[2]
    section = block * section_len
    ti1 = table_index - 1
    if verbose > 1:
        print("xi_tuple table", e1, ti1,  "; ", block, entry, source_shape)
        print("destination",  dest_tag , dest_shape) 
    ref_table = tables_xi.PERM_TABLES[e, table_index]
    ref_part = ref_table[section : section + section_len] 
    index = section + int(np.nonzero((ref_part & 0x7fff) == entry)[0][0])
    assert mm_sub_get_table_xi(e1, ti1, index, 0) & 0x7fff == entry
    if verbose > 1:
        print("dest tag", dest_tag, ", index", index)
    q, r = divmod(index, dest_shape[2])
    dest_index = 32 * q + r
    dest_sign = (mm_sub_get_table_xi(e1, ti1, q, 1) >> r) & 1
    ref_sign_table = tables_xi.SIGN_TABLES[e, table_index]
    assert dest_sign == (ref_sign_table[q] >> r) & 1
    return dest_sign, (dest_tag, dest_index)



def op_xi_table(v, exp, verbose):
    assert 1 <= exp <= 2
    exp1 = exp
    v_tuple = unit_vector_to_tuple(v)
    v_xi = tuple_to_xi_tuple(v_tuple)
    if verbose:
        print("xi tuple =", v_xi)
    sign, w_xi = op_xi_tuple(v_xi, exp1, verbose)
    w_tuple = xi_tuple_to_tuple(w_xi)
    if verbose:
       print("result tuple =", sign, w_tuple, ", xi_tuple =", w_xi)
    return (-1)**sign * space(*w_tuple)


def op_xi_ref(v, exp):
    g = group("l", exp)
    return v * g 


#####################################################################
# Tests
#####################################################################

def xi_testcases(n_cases = 50):
    test_cases = [
        (("B", 2, 1), 1),
        (("B", 12, 1), 1),
        (("C", 2, 1), 1),
        (("C", 12, 1), 1),
        (("T", 14, 1), 1),
        (("T", 21, 13), 2),
        (("T", 16, 1), 1),
        (("T", 18, 0), 1),
        (("T", 38, 10), 1),
        (("T", 376, 60), 1),
        (("T", 437, 63), 1),
        (("X", 0x400, 10), 2),
        (("X", 0x444, 10), 2),
    ]
    for v0, exp in test_cases:
        yield space(*v0), exp
    for i in range(n_cases // 2):
        i, j, e = randint(0,14), randint(0, 63), randint(1,2)
        yield space("T", i, j), e
    tags = "BCTX" 
    for v_tag in tags:  
        for i in range(n_cases):
            v = space(v_tag, "r")
            for exp in (1,2):
                yield v, exp
 

def one_test_op_xi(v, exp, verbose = 1):
    w_ref = op_xi_ref(v, exp)
    if verbose:
        print("\nCompute ", v, "* t**", exp, "=", w_ref)
    w = op_xi_table(v, exp, verbose)
    ok = w == w_ref
    if verbose or not ok:
        print(v, "* xi**", exp, "=")
        print(w)
        if not ok:
            print("expected:")
            print(w_ref)
            print("product with inverse:")
            print(op_xi_ref(v, 3-exp))
            raise ValueError("Monomial operation of xi failed")
  


@pytest.mark.mm_op
def test_op_xi(n_cases = 150, verbose = 0):
    print("Testing monomial operation of xi")
    i = 0
    for args in xi_testcases(n_cases):
        one_test_op_xi(*args, verbose = verbose)
        i += 1
    print("Test passed")





