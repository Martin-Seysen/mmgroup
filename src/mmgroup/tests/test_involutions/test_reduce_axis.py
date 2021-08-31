
from random import randint
import numpy as np

import pytest

from mmgroup import MM, MMSpace
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import leech2_matrix_basis
from mmgroup.clifford12 import leech2_matrix_radical
from mmgroup.clifford12 import leech2_matrix_expand
from mmgroup.clifford12 import leech_matrix_2A_axis_type
from mmgroup.mm import mm_aux_get_mmv1
from mmgroup.mm15 import op_word as mm_op15_word
from mmgroup.mm15 import op_find_short as mm_op15_find_short
from mmgroup.mm15 import op_t_A as mm_op15_t_A
from mmgroup.mm15 import op_compare as mm_op15_compare
from mmgroup.mm15 import op_store_axis as  mm_op15_store_axis
from mmgroup.mm15 import op_reduce_v_axis as  mm_op15_reduce_v_axis

V = MMSpace(15)


V_START_TUPLE = ("I", 3, 2)
V_START = V(V_START_TUPLE)
V_OPP = V_START * MM(('x', 0x200))

  
  
def leech_type(v2):
    """Return type of vector ``v2`` in the Leech lattice mod 2"""
    return gen_leech2_type(v2) >> 4




def short(v, value):
    """Return certain array of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the moster mod 15.

    The function returns the list of all vectors ``v2`` in the Leech
    lattice mod 2 such that the absolute value of the entry of ``v``
    corresponding a unit vector labelled by ``v2`` is equal to
    ``value``. That list is returned as a numpy array.
    """
    short = np.zeros(100000, dtype = np.uint32)
    l = mm_op15_find_short(v.data, short, len(short),  value, 0)  
    return short[:l]


def span(v, value, verbose = 1):
    """Return certain subspace of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the moster mod 15, and let
    ``L`` be the set of vectors in the Leech lattice mod 2
    returned by function ``short(v, value)``. The function returns
    linear subspace of the Leech lattice mod 2 spanned by ``L``
    as a list of vectors in a numpy array. 
    """
    short = np.zeros(100000, dtype = np.uint32)
    l = mm_op15_find_short(v.data, short, len(short),  value, 0)  
    short = short[:l]
    basis = np.zeros(24, dtype = np.uint64)
    if verbose:
        print("Function span, No of vectors =", l)
    dim= leech2_matrix_basis(short, l, basis, 24)
    if verbose:
        print("Dimension of spanned space:", dim)    
    basis = basis[:dim]
    v2new = np.zeros(1 << dim, dtype = np.uint32)
    l1 = leech2_matrix_expand(basis, dim, v2new)
    return v2new[:l1]

def radical(v, value, verbose = 1):
    """Return certain subspace of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the moster mod 15, and let
    ``L`` be the set of vectors in the Leech lattice mod 2
    returned by function ``short(v, value)``. let ``V`` be the
    linear space spanend by ``L``. 

    The function returns the radical of ``V`` (i.e. the intersection
    of ``V`` with its orthogonal complement) as a list of vectors in 
    a numpy array. 
    """
    short = np.zeros(100000, dtype = np.uint32)
    l = mm_op15_find_short(v.data, short, len(short),  value, 0)  
    short = short[:l]
    basis = np.zeros(24, dtype = np.uint64)
    if verbose:
        print("Function radical, No of vectors =", l)
    dim = leech2_matrix_radical(short, l, basis, 24)
    if verbose:
        print("Dimension of radical:", dim)    
    basis = basis[:dim]
    v2new = np.zeros(1 << dim, dtype = np.uint32)
    l1 = leech2_matrix_expand(basis, dim, v2new)
    return v2new[:l1]
       

##########################################################################


def reduce_axis(vector, verbose = 0):
    v = vector.data
    V = vector.space
    work = V()
    vA = np.zeros(24*4, dtype = np.uint64)
    r = np.zeros(120, dtype= np.uint32);
    len_r = 0
    for i in range(5):
        vt = leech_matrix_2A_axis_type(15, v)
        assert vt
        type = (vt >> 24) & 0xff
        vt &= 0xffffff
        if verbose: print("type",hex(type))
        if type == 0xC3:  # case 12C
            v2all = radical(v, 7, verbose)  
            v4all = [w for w in v2all if leech_type(w) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x42, 0x61]
        elif type == 0xA2:  # case 10B
            v2all = radical(v, 4, verbose)
            v4all = [w for w in v2all if leech_type(w) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x42, 0x43]
        elif type == 0xA1:  # case 10A
            v0 = short(v, 3)[0]
            v2all = short(v, 1)
            v4all = [w ^ v0 for w in v2all if leech_type(w ^ v0) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x61]
        elif type == 0x82:  # case 8B
            v2all = short(v, 1)
            v2_0 = v2all[0]
            v4all = [w ^ v2_0 for w in v2all if leech_type(w ^ v2_0) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x41]
        elif type == 0x66:  # case 6F
            v2all = radical(v, 7, verbose)
            v4all = [w for w in v2all if leech_type(w) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x43]
        elif type == 0x63:  # case 6C
            v2all = span(v, 3, verbose)
            v4all = [w  for w in v2all if leech_type(w) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x41]
        elif type == 0x61:  # case 6A
            v2all = short(v, 5)
            v4all = [w ^ vt for w in v2all if leech_type(w ^vt) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x41]
        elif type == 0x43:  # case 4C
            v2all = radical(v, 1, verbose)
            v4all = [w for w in v2all if leech_type(w) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x22]
        elif type == 0x42:  # case 4B
            v2all = radical(v, 1, verbose)
            v4all = [w for w in v2all if leech_type(w) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x22]
        elif type == 0x41:  # case 4A
            v4 = vt
            t_types = [0x21]
        elif type == 0x22:  # case 2B
            v2all = span(v, 4, verbose)  
            v4all = [w for w in v2all if leech_type(w) == 4]
            assert len(v4all)
            v4 = v4all[0]
            t_types = [0x21]
        elif type == 0x21:  # case 2A
            r1 = gen_leech2_reduce_type2(vt, 0, r[len_r:])
            assert r1 >= 0
            mm_op15_word(v, r[len_r:], r1, 1, work.data)
            len_r += r1
            ind = mm_aux_get_mmv1(15, v, (24+3)*32 + 2)
            if (ind != 15-2):
                r[len_r] = 0xB0000200
                mm_op15_word(v, r[len_r:], 1, 1, work.data)
                len_r += 1
            assert mm_op15_compare(v, V_START.data) == 0
            return r[:len_r]
        else:
            raise ValueError("WTF1")

        r1 = gen_leech2_reduce_type4(v4, r[len_r:])
        assert r1 >= 0
        mm_op15_word(v, r[len_r:], r1, 1, work.data)
        len_r += r1
        ok = False
        for e in (1,2):
            mm_op15_t_A(v, e, vA)                    
            t = leech_matrix_2A_axis_type(15, vA) >> 24
            if verbose: print("e", e, hex(t))
            if t in t_types:
                r[len_r] = 0xD0000003 - e
                mm_op15_word(v, r[len_r:], 1, 1, work.data)
                len_r += 1
                ok = True
                break
        if not ok:
            raise ValueError("WTF2")

    raise ValueError("WTF3")    
        
        
    







##############################################################################


from mmgroup.tests.test_involutions.test_2A_axes import AXES




def make_testcases():
    V = V_START.space
    yield V_START.copy()
    yield V(  ("I", 11, 9) )
    for i in range(10):
        yield V_START * MM.rand_G_x0()
    for ax in AXES:
        #if ax in ["10B", "12C"]: continue   ## preliminary!!!!!!!
        v0 = V_START * MM(AXES[ax])
        for i in range(20):
            yield v0 * MM.rand_G_x0()
    for quality in range(2,11):
        for i in range(3):
              yield  V_START *  MM.rand_mm(quality)      


@pytest.mark.involution
def test_reduce_axis(verbose = 0):
    for i, v in enumerate(make_testcases()):
        if verbose:
            print("Test case", i)
        r = reduce_axis(v.copy())
        g = MM.from_data(r)
        assert v * g == V_START

        vr1 = np.zeros(200, dtype = np.uint32)
        len_r1 = mm_op15_reduce_v_axis(v.copy().data, vr1)
        assert len_r1 >= 0
        g1 = MM.from_data(vr1[:len_r1])
        assert g1 == g





