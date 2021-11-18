
from random import randint
import numpy as np

import datetime
import time
import pytest

from mmgroup import MM0, MMSpace, MMV
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type2_ortho
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import leech2_matrix_basis
from mmgroup.clifford12 import leech2_matrix_radical
from mmgroup.clifford12 import leech2_matrix_expand
from mmgroup.mm import mm_aux_get_mmv1
from mmgroup.mm15 import op_word as mm_op15_word
from mmgroup.mm15 import op_axes_find_short as mm_op15_axes_find_short
from mmgroup.mm15 import op_t_A as mm_op15_t_A
from mmgroup.mm15 import op_compare as mm_op15_compare
from mmgroup.mm15 import op_store_axis as  mm_op15_store_axis
from mmgroup.mm15 import op_reduce_v_axis as  mm_op15_reduce_v_axis
from mmgroup.mm15 import op_reduce_v_baby_axis as mm_op15_reduce_v_baby_axis
from mmgroup.mm15 import op_reduce_G_x0 as mm_op15_reduce_G_x0
from mmgroup.mm15 import op_2A_axis_type as mm_op15_2A_axis_type
from mmgroup.mm15 import op_eval_A as mm_op15_eval_A

from mmgroup.tests.test_axes.test_reduce_type2 import rand_Co2
from mmgroup.tests.test_axes.test_import import AXES, BABY_AXES

V = MMV(15)


V_START_TUPLE = ("I", 3, 2)
V_START = V(*V_START_TUPLE)
V_OPP = V_START * MM0('x', 0x200)

v_start = 0x200


#########################################################################
# Auxiliary functions 
#########################################################################
 
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
    l = mm_op15_axes_find_short(v.data, short, len(short),  value, 0)  
    return short[:l]


def span(v, value, verbose = 0):
    """Return certain subspace of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the moster mod 15, and let
    ``L`` be the set of vectors in the Leech lattice mod 2
    returned by function ``short(v, value)``. The function returns
    linear subspace of the Leech lattice mod 2 spanned by ``L``
    as a list of vectors in a numpy array. 
    """
    short = np.zeros(100000, dtype = np.uint32)
    l = mm_op15_axes_find_short(v.data, short, len(short),  value, 0)  
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

def radical(v, value, verbose = 0):
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
    l = mm_op15_axes_find_short(v.data, short, len(short),  value, 0)  
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
       


def find_type4(v_list):
    for w in v_list:
        if leech_type(w) == 4:
            return w
    err = "No type-4 vector found"
    raise ValueError(err)
  

def find_ortho_short(vlist):
    for w in vlist:
        if leech_type(w) == 4 and leech_type(w ^ v_start) == 2:
            return w ^ v_start
    err = "No short vector orthogonal to 'v_start' found"
    raise ValueError(err)
  
def eval_A_vstart(v):
    return mm_op15_eval_A(v, 0x200)

##########################################################################
# Reducing a 2A axis to V_START
##########################################################################


def reduce_axis(vector, verbose = 0):
    v = vector.data
    #V = vector.space
    work = V()
    vA = np.zeros(24*4, dtype = np.uint64)
    r = np.zeros(120, dtype= np.uint32);
    len_r = 0
    if verbose: print("Function reduce_axis")
    for i in range(5):
        vt = mm_op15_2A_axis_type(v)
        assert vt
        type = (vt >> 24) & 0xff
        vt &= 0xffffff
        if verbose: print("type",hex(type))
        if type == 0xC3:  # case 12C
            v2all = radical(v, 7, verbose)  
            v4 = find_type4(v2all)
            t_types = [0x42, 0x61]
        elif type == 0xA2:  # case 10B
            v2all = radical(v, 4, verbose)
            v4 = find_type4(v2all)
            t_types = [0x42, 0x43]
        elif type == 0xA1:  # case 10A
            v0 = short(v, 3)[0]
            v2all = short(v, 1)
            v2all = [w ^ v0 for w in v2all]
            v4 = find_type4(v2all)
            t_types = [0x61]
        elif type == 0x82:  # case 8B
            v2all = short(v, 1)
            v2_0 = v2all[0]
            v2all = [w ^ v2_0 for w in v2all] 
            v4 = find_type4(v2all)
            t_types = [0x41]
        elif type == 0x66:  # case 6F
            v2all = radical(v, 7, verbose)
            v4 = find_type4(v2all)
            t_types = [0x43]
        elif type == 0x63:  # case 6C
            v2all = span(v, 3, verbose)
            v4 = find_type4(v2all)
            t_types = [0x41]
        elif type == 0x61:  # case 6A
            v2all = short(v, 5)
            v2all = [w ^ vt for w in v2all]
            v4 = find_type4(v2all)
            t_types = [0x41]
        elif type == 0x43:  # case 4C
            v2all = radical(v, 1, verbose)
            v4 = find_type4(v2all)
            t_types = [0x22]
        elif type == 0x42:  # case 4B
            v2all = radical(v, 1, verbose)
            v4 = find_type4(v2all)
            t_types = [0x22]
        elif type == 0x41:  # case 4A
            v4 = vt
            t_types = [0x21]
        elif type == 0x22:  # case 2B
            v2all = span(v, 4, verbose)  
            v4 = find_type4(v2all)
            t_types = [0x21]
        elif type == 0x21:  # case 2A
            r1 = gen_leech2_reduce_type2(vt, r[len_r:])
            assert r1 >= 0
            mm_op15_word(v, r[len_r:], r1, 1, work.data)
            len_r += r1
            ind = mm_aux_get_mmv1(15, v, (24+3)*32 + 2)
            if (ind != 15-2):
                r[len_r] = 0xB0000200
                mm_op15_word(v, r[len_r:], 1, 1, work.data)
                len_r += 1
            assert mm_op15_compare(v, V_START.data) == 0
            if verbose: 
                print("Function reduce_axis terminated successfullly")
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
            t = mm_op15_2A_axis_type(vA) >> 24
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
        
        
    

##########################################################################
# Reducing a 2A axis to orthogonal to V_START, preserving V_START
##########################################################################



def reduce_baby_axis(vector, verbose = 1):
    v = vector.data
    #V = vector.space
    work = V()
    vA = np.zeros(24*4, dtype = np.uint64)
    r = np.zeros(120, dtype= np.uint32);
    len_r = 0
    if verbose: print("Function reduce_baby_axis")
    for i in range(5):
        vt = mm_op15_2A_axis_type(v)
        assert vt
        type = (vt >> 24) & 0xff
        vt &= 0xffffff
        if verbose:
            print("type =", hex(type), ", vt =", hex(vt), 
                 ", A(v0) =",  eval_A_vstart(v))
        if type == 0xA1:  # case 10A
            v0 = short(v, 3)[0]
            v2all = short(v, 1)
            v2all = [w ^ v0 for w in v2all]
            v4 = find_ortho_short(v2all)
            t_types = [0x61]
        elif type == 0x63:  # case 6C
            v2all = span(v, 3, verbose)
            v4 = find_ortho_short(v2all)
            t_types = [0x41]
        elif type == 0x61:  # case 6A
            v2all = short(v, 5)
            v2all = [w ^ vt for w in v2all] # if leech_type(w ^vt) == 4]
            v4 = find_ortho_short(v2all)
            t_types = [0x41, 0x43]
        elif type == 0x43:  # case 4C
            v2all = radical(v, 1, verbose)
            v4 = find_ortho_short(v2all)
            t_types = [0x22]
        elif type == 0x42:  # case 4B
            v2all = radical(v, 1, verbose)
            v4 = find_ortho_short(v2all)
            t_types = [0x22]
        elif type == 0x41:  # case 4A
            v4 = vt ^ v_start
            assert leech_type(v4) == 2
            t_types = [0x21]
        elif type == 0x22:  # case 2B
            if eval_A_vstart(v) in [0, 8]:
                v2all = span(v, 4, verbose)  
                v4 = find_ortho_short(v2all)
                t_types = [0x21]
            else:
                raise ValueError("WTF 2B")
        elif type == 0x21:  # case 2A
            if eval_A_vstart(v) == 0:
                r1 = gen_leech2_reduce_type2_ortho(vt, r[len_r:])
                assert r1 >= 0
                mm_op15_word(v, r[len_r:], r1, 1, work.data)
                len_r += r1
                vt = mm_op15_2A_axis_type(v) & 0xffffff
                assert vt == 0x800200 
                ind = mm_aux_get_mmv1(15, v, (2*24+3)*32 + 2)
                e = 2 - (ind == 15-2)
                if verbose: print("ind", ind, e)
                r[len_r] = 0xD0000003 - e
                mm_op15_word(v, r[len_r:], 1, 1, work.data)
                len_r += 1
            assert  eval_A_vstart(v) == 4, eval_A_vstart(v)
            if verbose: 
                print("Function reduce_baby_axis terminated successfullly")
            return r[:len_r]
        else:
            raise ValueError("WTF1")

        r1 = gen_leech2_reduce_type2_ortho(v4, r[len_r:])
        assert r1 >= 0
        mm_op15_word(v, r[len_r:], r1, 1, work.data)
        len_r += r1
        ok = False
        for e in (1,2):
            mm_op15_t_A(v, e, vA)                    
            t = mm_op15_2A_axis_type(vA) >> 24
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
        
        
    




##########################################################################
# Testing function reduce_axis and C function mm_op15_reduce_v_axis
##########################################################################







def make_axis_testcases():
    #V = V_START.space
    yield V_START.copy()
    yield V("I", 11, 9)
    for i in range(10):
        yield V_START * MM0("r", "G_x0")
    for ax in AXES:
        v0 = V_START * MM0(AXES[ax])
        for i in range(20):
            yield v0 * MM0("r", "G_x0")
    for quality in range(2,11):
        for i in range(3):
              yield  V_START *  MM0("r", quality)      



@pytest.mark.axes
def test_reduce_axis(verbose = 0):
    for i, v in enumerate(make_axis_testcases()):
        if verbose:
            print("\nTest case", i)
        r = reduce_axis(v.copy(), verbose)
        g = MM0('a', r)
        assert v * g == V_START

        vr1 = np.zeros(200, dtype = np.uint32)
        len_r1 = mm_op15_reduce_v_axis(v.copy().data, vr1)
        assert len_r1 >= 0
        g1 = MM0('a', vr1[:len_r1])
        assert g1 == g




##########################################################################
# Testing function reduce_baby_axis and mm_op15_reduce_v_baby_axis
##########################################################################

def rand_BM(quality = 8):
    r"""Generate certain 'random' element in the Baby Monster`

    Let ``v_0`` be the element of the subgroup  ``Q_x0`` of ``G_x0``
    corresponding to the Golay cocode element ``[2,3]``.

    The function generates a random element of the centralizer
    of ``v_0`` in the monster.
    """
    a = MM0()
    for i in range(quality):
         e = randint(0, 2)
         a *= MM0([rand_Co2(), ('t', e)]) 
    return a 




def make_baby_testcases():
    #V = V_START.space
    yield V_OPP.copy()
    for i in range(10):
        yield V_OPP.copy() * rand_Co2()
    for ax in BABY_AXES:
        v0 = V_OPP * MM0(BABY_AXES[ax])
        for i in range(5):
            yield v0 * rand_BM()
    for quality in range(2,11):
        for i in range(5):
              yield  V_OPP.copy() * rand_BM(quality)      



@pytest.mark.axes
def test_reduce_baby_axis(verbose = 0):
    for i, v in enumerate(make_baby_testcases()):
        if verbose:
            print("\nTest case", i)
        r = reduce_baby_axis(v.copy(), verbose)
        g = MM0('a', r)
        assert v * g == V_OPP
        assert V_START * g == V_START

        vr1 = np.zeros(200, dtype = np.uint32)

        len_r1 = mm_op15_reduce_v_baby_axis(v.copy().data, vr1)
        assert len_r1 >= 0
        g1 = MM0('a', vr1[:len_r1])

        ok =  g1 == g 
        #print(type(g1), type(g))
        if verbose or not ok:
             vt = mm_op15_2A_axis_type(v.data) >> 24
             vA = eval_A_vstart(v.data)
             print("Type(v) = 0x%x, value(A) = %d" % (vt, vA))
             print("Op:  ", [hex(x) for x in g] )
             print("Fast:", [hex(x) for x in g1] ) 
             if not ok:
                 err = "Function mm_op15_reduce_v_baby_axis failed"
                 raise ValueError(err)


##########################################################################
# Testing C function mm_op15_reduce_G_x0
##########################################################################


reduce_time = None

def reduce_G_x0(g):
    global reduce_time
    #assert g.group == MM
    r = np.zeros(256, dtype = np.uint32)
    t_start = time.perf_counter() 
    res = mm_op15_reduce_G_x0(g.mmdata, len(g.mmdata), r)
    reduce_time = time.perf_counter() - t_start
    if res < 0:
       err = "Error %d in reducing element of monster group"
       raise ValueError(err % res)
    return MM0('a', r[:res])


def g_complexity(g):
    s = [x for x in g.mmdata if x & 0x70000000 == 0x50000000]
    return len(s)


def make_reduce_testcases():
    for quality in range(1,17):
        for i in range(2):
              yield  MM0('r', quality)      


@pytest.mark.axes
def test_reduce_G_x0(verbose = 0):
    for i, g in enumerate(make_reduce_testcases()):
        if verbose:
            print("\nTest case", i)
        g1 = reduce_G_x0(g)
        g2 = g * g1
        ok = g2.in_G_x0()

        if verbose or not ok:
            print("g =", g)
            print("g1 =", g1 ) 
            print("g2 =", g2 ) 
            print("Run time: %.2f ms, Complexity: %d" % 
                (1000 * reduce_time, g_complexity(g))) 
            if not ok:
                err = "Function mm_op15_reduce_G_x0 failed"
                raise ValueError(err)



