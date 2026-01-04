from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from random import randint

import pytest


from mmgroup.mm_op import mm_sub_get_table_xi
from mmgroup.mm_op import mm_sub_get_offset_table_xi
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



MAP_OFFSET_BOX = {
   0x300 : 'BC', 0xcc0 : 'T0', 0x66c0 : 'T1', 0xc6c0 : 'X0', 0x146c0: 'X1'
}


def source_box_to_table_index(box, e):
    e1 = (e - 1) % 3
    assert e1 < 2
    for i in range(5):
        ofs = mm_sub_get_offset_table_xi(i, e1, 0)
        #print(box, e, i, e1, 0, hex(ofs))
        if ofs in MAP_OFFSET_BOX and MAP_OFFSET_BOX[ofs] == box:
            return i
    ERR = "Could not find table index for box %s, exponent %d"
    raise ValueError(ERR % (box, e))

def table_index_to_boxes(index, e):
    e1 = (e - 1) % 3
    assert e1 < 2
    assert 0 <= index < 5
    src_box = MAP_OFFSET_BOX[mm_sub_get_offset_table_xi(index, e1, 0)]
    dest_box = MAP_OFFSET_BOX[mm_sub_get_offset_table_xi(index, e1, 1)]
    return src_box, dest_box


def op_xi_tuple(v, e, verbose =  0):
    from mmgroup.dev.mm_basics.mm_tables_xi import MM_TablesXi
    tables_xi = MM_TablesXi()
    e1 = e - 1
    """
    table_index = None
    for i, source_tag in tables_xi.SOURCE_TAGS.items():
        if source_tag[e1] == v[0]:
            table_index = i
    dest_tag = tables_xi.DEST_TAGS[table_index][e1]
    """
    index = source_box_to_table_index(v[0], e)
    #print("index =", index)
    source_box, dest_box = table_index_to_boxes(index, e)
    assert source_box == v[0]
    source_shape, dest_shape = tables_xi.SHAPES[index]
    block, entry = divmod(v[1], source_shape[1] * 32)
    section_len = dest_shape[1] * dest_shape[2]
    section = block * section_len
    if verbose > 1:
        print("xi_tuple table", e1, index,  "; ", block, entry, source_shape)
        print("destination",  dest_box , dest_shape) 
    ref_table = tables_xi.PERM_TABLES[index][e1]
    ref_part = ref_table[section : section + section_len] 
    j = section + int(np.nonzero((ref_part & 0x7fff) == entry)[0][0])
    assert mm_sub_get_table_xi(index, e1, j, 0) & 0x7fff == entry
    if verbose > 1:
        print("dest tag", dest_box, ", index", j)
    q, r = divmod(j, dest_shape[2])
    dest_index = 32 * q + r
    dest_sign = (mm_sub_get_table_xi(index, e1, q, 1) >> r) & 1
    ref_sign_table = tables_xi.SIGN_TABLES[index][e1]
    assert dest_sign == (ref_sign_table[q] >> r) & 1
    return dest_sign, (dest_box, dest_index)



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
 

def one_test_op_xi(v, exp, verbose = 0):
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



#####################################################################
# In the remainder of this file we test the following program
#####################################################################






def map_xi(v_in, e, v_out):
        r"""Compute v_out = v_in * \xi**e, for e = 1, 2.

        We compute the monomial part of v_out only.
        v_in is a array of integers corresponding to a vector in the
        representation \rho of the Monster. Vector v_out is
        an empty vector of integers of the same size.
        """
        for stage in range(5):
            box_in = v_in[OFFSETS[stage][e-1][0]:]
            box_out = v_out[OFFSETS[stage][e-1][1]:]
            shape_in = SHAPES[stage][0]
            shape_out = SHAPES[stage][1]
            cluster_in_size = shape_in[1] * 32
            cluster_out_size = shape_out[1] * 32
            cluster_perm_size = shape_out[1] * shape_out[2]
            cluster_sign_size = shape_out[1]
            perm_table = PERM_TABLES[stage][e-1]
            sign_table = SIGN_TABLES[stage][e-1]
            for cluster in range(shape_out[0]):
                cluster_in = box_in[cluster * cluster_in_size:]
                cluster_out = box_out[cluster * cluster_out_size:]
                cluster_perm = perm_table[cluster * cluster_perm_size:]
                cluster_sign = sign_table[cluster * cluster_sign_size:]
                for i in range(shape_out[1]):
                    for j in range(shape_out[2]):
                        x = int(cluster_in[cluster_perm[shape_out[2]*i + j]])
                        x = (-1)**(int(cluster_sign[i]) >> j) * x
                        cluster_out[32 * i + j] = x



from mmgroup.mm_op import mm_aux_index_extern_to_intern
from mmgroup.mm_op import mm_aux_index_intern_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
from mmgroup.mm_op import mm_aux_index_intern_to_leech2
from mmgroup.mm_op import mm_aux_index_leech2_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_intern
from mmgroup.mm_op import mm_aux_index_leech2_to_intern_fast
from mmgroup.generators import gen_xi_op_xi


def load_tables():
    from mmgroup.dev.mm_basics.mm_tables_xi import MM_TablesXi
    global PERM_TABLES, SIGN_TABLES, OFFSETS, SHAPES, LEN_VECTOR
    tbl = MM_TablesXi()
    PERM_TABLES = tbl.PERM_TABLES
    SIGN_TABLES = tbl.SIGN_TABLES
    OFFSETS = tbl.OFFSETS
    SHAPES = tbl.SHAPES
    LEN_VECTOR = 247488



def demo_xi_test_indices():
    for i in range(300, 300+98280, 111):
         yield  mm_aux_index_extern_to_intern(i)


def map_xi_intern(x, e):
    x1 = mm_aux_index_sparse_to_leech2(mm_aux_index_intern_to_sparse(x))
    assert x1 == mm_aux_index_intern_to_leech2(x)
    x2 = mm_aux_index_leech2_to_intern_fast(x1)
    assert x == x2, (hex(x), hex(x1), hex(x2))
    # print(hex(x), hex(x1))
    y = gen_xi_op_xi(x1 & 0xffffff, e)
    sign = (-1) ** (y >> 24)
    y = mm_aux_index_sparse_to_intern(mm_aux_index_leech2_to_sparse(y))
    return y, sign



@pytest.mark.mm_op
def test_demo_xi_monomimal():
    load_tables()
    for e in (1,2):
        v_in = np.random.randint(-127, 127, LEN_VECTOR, dtype=np.int8)
        v_out = np.zeros(LEN_VECTOR, dtype=np.int8)
        map_xi(v_in, e, v_out)
        for index_in in demo_xi_test_indices():
             entry_in = v_in[index_in]
             index_out, sign = map_xi_intern(index_in, e)
             entry_out = v_out[index_out]
             if entry_out !=  sign * entry_in:
                 print("\nMap index %s with eponent %d" % (hex(index_in), e))
                 print("Entry is", entry_in)
                 print("Output index = %s, sign = %d" % (hex(index_out), sign))
                 print("Output entry obtained:", entry_out)
                 ERR = "Test'test_demo_xi_monomimal' failed"
                 raise ValueError(ERR)
                 #print(ERR)

