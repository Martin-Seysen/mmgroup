
from random import randint, shuffle
import numpy as np

import datetime
import time
import pytest

from mmgroup import mat24
from mmgroup import MM0, MMSpace, MMV, Cocode, AutPL, XLeech2
from mmgroup.generators import gen_leech2_mul
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type2_ortho
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_op_word
from mmgroup.clifford12 import leech2_matrix_basis
from mmgroup.clifford12 import leech2_matrix_radical
from mmgroup.clifford12 import leech2_matrix_expand
from mmgroup.mm_op import mm_vector
from mmgroup.mm_op import mm_aux_get_mmv, mm_aux_get_mmv_leech2
from mmgroup.mm_op import mm_op_word
from mmgroup.mm_op import mm_op_eval_X_find_abs
from mmgroup.mm_op import mm_op_t_A 
from mmgroup.mm_op import mm_op_compare
from mmgroup.mm_op import mm_op_copy
from mmgroup.mm_op import mm_op_eval_A


def import_all():
    global mm_reduce_find_type4
    global mm_reduce_2A_axis_type
    global mm_reduce_load_axis
    global mm_reduce_vector_vp
    global mm_reduce_vector_vm
    global mm_reduce_vector_shortcut
    global get_baby_sample_axes
    global get_sample_axes

    from mmgroup.mm_reduce import mm_reduce_find_type4
    from mmgroup.mm_reduce import mm_reduce_2A_axis_type
    from mmgroup.mm_reduce import mm_reduce_load_axis
    from mmgroup.mm_reduce import mm_reduce_vector_vp
    from mmgroup.mm_reduce import mm_reduce_vector_vm
    from mmgroup.mm_reduce import mm_reduce_vector_shortcut
    from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
    from mmgroup.tests.axes.get_sample_axes import get_sample_axes

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
    return gen_leech2_type(v2)


S_AXIS = "  %d coordinates of the axis in 98280_x have absolute value %d"

def short(v, value, verbose = 0):
    """Return certain array of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the monster mod 15.

    The function returns the list of all vectors ``v2`` in the Leech
    lattice mod 2 such that the absolute value of the entry of ``v``
    corresponding a unit vector labelled by ``v2`` is equal to
    ``value``. That list is returned as a numpy array.
    """
    short = np.zeros(100000, dtype = np.uint32)
    l = mm_op_eval_X_find_abs(15, v.data, short, len(short),  value, 0) 
    if verbose:
        print(S_AXIS % (l, value))
    return short[:l]


def span(v, value, verbose = 0):
    """Return certain subspace of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the monster mod 15, and let
    ``L`` be the set of vectors in the Leech lattice mod 2
    returned by function ``short(v, value)``. The function returns
    linear subspace of the Leech lattice mod 2 spanned by ``L``
    as a list of vectors in a numpy array. 
    """
    short = np.zeros(100000, dtype = np.uint32)
    l = mm_op_eval_X_find_abs(15, v.data, short, len(short), value, 0)  
    short = short[:l]
    basis = np.zeros(24, dtype = np.uint64)
    dim = leech2_matrix_basis(short, l, basis, 24)
    if verbose:
        #print(S_AXIS % (l, value))
        print("  Dimension of spanned space:", dim)    
    basis = basis[:dim]
    v2new = np.zeros(1 << dim, dtype = np.uint32)
    l1 = leech2_matrix_expand(basis, dim, v2new)
    return v2new[:l1]

def radical(v, value, verbose = 0):
    """Return certain subspace of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the monster mod 15, and let
    ``L`` be the set of vectors in the Leech lattice mod 2
    returned by function ``short(v, value)``. let ``V`` be the
    linear space spanend by ``L``. 

    The function returns the radical of ``V`` (i.e. the intersection
    of ``V`` with its orthogonal complement) as a list of vectors in 
    a numpy array. 
    """
    short = np.zeros(100000, dtype = np.uint32)
    l = mm_op_eval_X_find_abs(15, v.data, short, len(short),  value, 0)  
    short = short[:l]
    basis = np.zeros(24, dtype = np.uint64)
    dim = leech2_matrix_radical(short, l, basis, 24)
    if verbose:
        print(S_AXIS % (l, value))
        print("  Dimension of radical:", dim)    
    basis = basis[:dim]
    v2new = np.zeros(1 << dim, dtype = np.uint32)
    l1 = leech2_matrix_expand(basis, dim, v2new)
    return v2new[:l1]
       


def find_type4_old(v_list):
    for w in v_list:
        if leech_type(w) == 4:
            return w
    err = "No type-4 vector found"
    raise ValueError(err)
  

def find_type4(v_list):
    a = np.array(v_list, dtype = np.uint32)
    v = mm_reduce_find_type4(a, len(a), 0)
    if v == 0:
        raise ValueError("No type-4 vector found in list")
    return v

def find_ortho_short_old(vlist):
    for w in vlist:
        if leech_type(w) == 4 and leech_type(w ^ v_start) == 2:
            return w ^ v_start
    err = "No short vector orthogonal to 'v_start' found"
    raise ValueError(err)
  
def find_ortho_short(v_list, axis = v_start):
    a = np.array(v_list, dtype = np.uint32)
    v = mm_reduce_find_type4(a, len(a), axis)
    if v == 0:
        raise ValueError("No feasible type-4 vector found in list")
    return v





def v_leech2_adjust_sign(v, v2):
    v2 &= 0xffffff
    value = mm_aux_get_mmv_leech2(15, v, v2)
    assert value in (2, 13)
    if (value == 2): v2 += 0x1000000
    return v2


def eval_A_vstart(v, v_start = 0x200):
    return mm_op_eval_A(15, v, v_start)

##########################################################################
# Reducing a 2A axis to V_START
##########################################################################


def reduce_axis(vector, std_axis = 1, verbose = 0):
    v = vector.data
    #V = vector.space
    work = V()
    vA = np.zeros(24*4, dtype = np.uint64)
    r = np.zeros(120, dtype= np.uint32);
    len_r = 0
    if verbose: print("Function reduce_axis")
    for i in range(5):
        vt = mm_reduce_2A_axis_type(v)
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
            if not std_axis:
                vt = v_leech2_adjust_sign(v, vt)
                if verbose: 
                    print("Function reduce_axis terminated successfullly")
                return r[:len_r], vt
            r1 = gen_leech2_reduce_type2(vt, r[len_r:])
            assert r1 >= 0
            mm_op_word(15, v, r[len_r:], r1, 1, work.data)
            len_r += r1
            ind = mm_aux_get_mmv(15, v, (24+3)*32 + 2)
            if (ind != 15-2):
                r[len_r] = 0xB0000200
                mm_op_word(15, v, r[len_r:], 1, 1, work.data)
                len_r += 1
            assert mm_op_compare(15, v, V_START.data) == 0
            if verbose: 
                print("Function reduce_axis terminated successfullly")
            return r[:len_r], 0x200
        else:
            raise ValueError("WTF1")

        r1 = gen_leech2_reduce_type4(v4, r[len_r:])
        assert r1 >= 0
        mm_op_word(15, v, r[len_r:], r1, 1, work.data)
        len_r += r1
        ok = False
        for e in (1,2):
            mm_op_t_A(15, v, e, vA)                    
            t = mm_reduce_2A_axis_type(vA) >> 24
            if verbose: print("e", e, hex(t))
            if t in t_types:
                r[len_r] = 0xD0000003 - e
                mm_op_word(15, v, r[len_r:], 1, 1, work.data)
                len_r += 1
                ok = True
                break
        if not ok:
            raise ValueError("WTF2")

    raise ValueError("WTF3")    
        
        
    

##########################################################################
# Reducing a 2A axis to orthogonal to V_START, preserving V_START
##########################################################################



def op_axis(axis, r):
    transformed_axis = gen_leech2_op_word(axis, r, len(r))
    if axis == v_start:
        assert transformed_axis == axis;
    return transformed_axis; 


def reduce_baby_axis(vector, axis = v_start, verbose = 0):
    v = vector.data
    #V = vector.space
    work = V()
    vA = np.zeros(24*4, dtype = np.uint64)
    r = np.zeros(120, dtype= np.uint32)
    len_r = 0
    if verbose: print("Function reduce_baby_axis")
    for i in range(5):
        vt = mm_reduce_2A_axis_type(v)
        assert vt
        type = (vt >> 24) & 0xff
        vt &= 0xffffff
        if verbose:
            print("type =", hex(type), ", vt =", hex(vt), 
                 ", A(v0) =",  mm_op_eval_A(15, v, axis), 
                 ", axis =", hex(axis))
        if type == 0xA1:  # case 10A
            v0 = short(v, 3)[0]
            v2all = short(v, 1)
            v2all = [w ^ v0 for w in v2all]
            v4 = find_ortho_short(v2all, axis)
            t_types = [0x61]
        elif type == 0x63:  # case 6C
            v2all = span(v, 3, verbose)
            v4 = find_ortho_short(v2all, axis)
            t_types = [0x41]
        elif type == 0x61:  # case 6A
            v2all = short(v, 5)
            v2all = [w ^ vt for w in v2all] # if leech_type(w ^vt) == 4]
            v4 = find_ortho_short(v2all, axis)
            t_types = [0x41, 0x43]
        elif type == 0x43:  # case 4C
            v2all = radical(v, 1, verbose)
            v4 = find_ortho_short(v2all, axis)
            t_types = [0x22]
        elif type == 0x42:  # case 4B
            v2all = radical(v, 1, verbose)
            v4 = find_ortho_short(v2all, axis)
            t_types = [0x22]
        elif type == 0x41:  # case 4A
            v4 = vt
            assert leech_type(v4) == 4
            t_types = [0x21]
        elif type == 0x22:  # case 2B
            if mm_op_eval_A(15, v, axis) in [0, 8]:
                v2all = span(v, 4, verbose)  
                v4 = find_ortho_short(v2all, axis)
                t_types = [0x21]
            else:
                raise ValueError("WTF 2B")
        elif type == 0x21:  # case 2A
            if verbose: 
                print("v has type 0x21, axis = %s, eval_A = %d" %
                     (hex(axis),  mm_op_eval_A(15, v, axis)) )
            if mm_op_eval_A(15, v, axis) == 0:
                vt ^= axis
                r1 = gen_leech2_reduce_type4(vt, r[len_r:])
                assert r1 >= 0
                mm_op_word(15, v, r[len_r:], r1, 1, work.data)
                axis = op_axis(axis,  r[len_r: len_r + r1])
                len_r += r1
                vt = mm_reduce_2A_axis_type(v) & 0xffffff
                vt = v_leech2_adjust_sign(v, vt)
                vt = gen_leech2_mul(vt, axis)
                assert vt & 0xffffff == 0x800000 
                e = 1 + (vt >> 24);
                if verbose: print("e",  e)
                r[len_r] = 0xD0000003 - e
                mm_op_word(15, v, r[len_r:], 1, 1, work.data)
                axis = op_axis(axis,  r[len_r: len_r + 1])
                len_r += 1
            else:
                assert mm_op_eval_A(15, v, axis) == 4
            vt = mm_reduce_2A_axis_type(v) & 0xffffff
            vt = v_leech2_adjust_sign(v, vt)
            assert (vt ^ axis) & 0x1ffffff == 0x1000000, (hex(vt), hex(axis))
            if verbose: 
                print("Function reduce_baby_axis terminated successfullly")
            return r[:len_r], axis
        else:
            raise ValueError("WTF1")

        r1 = gen_leech2_reduce_type4(v4, r[len_r:])
        assert r1 >= 0, (hex(v4), )
        mm_op_word(15, v, r[len_r:], r1, 1, work.data)
        axis = op_axis(axis,  r[len_r: len_r + r1])
        len_r += r1
        ok = False
        for e in (1,2):
            mm_op_t_A(15, v, e, vA)                    
            t = mm_reduce_2A_axis_type(vA) >> 24
            if verbose: print("e", e, hex(t))
            if t in t_types:
                r[len_r] = 0xD0000003 - e
                mm_op_word(15, v, r[len_r:], 1, 1, work.data)
                axis = op_axis(axis,  r[len_r: len_r + 1])
                len_r += 1
                ok = True
                break
        if not ok:
            raise ValueError("WTF2")

    raise ValueError("WTF3")    
        
        
    

#########################################################################
## Generating a random vector in the subgroup ``Co_2`` of ``G_x0``
#########################################################################

v_start_0 =  0x200 
assert Cocode(v_start_0).syndrome_list() == [2,3]
v_ortho_start = 0x800200


def rand_pi_mat22():
    r"""Generate certain 'random' element of AutPL

    The function generates a random element ``e`` of the automorphism
    group ``AutPL`` of the Parker loop, such that the permutation
    in the Mathieu ``M24`` group corresponding to ``e`` preserves the 
    set  ``\{2, 3\}``.
    """
    pi = AutPL('r', 'r')
    perm = pi.perm
    l_src = perm[2:4]
    shuffle(l_src)
    _, a = mat24.perm_from_map(l_src, [2, 3])
    return  pi  *  AutPL('r', a)  

def rand_Co2(quality = 5):
    r"""Generate certain 'random' element in a subgroup of ``G_x0``

    Let ``v_0`` be the element of the subgroup  ``Q_x0`` of ``G_x0``
    corresponding to the Golay cocode element ``[2,3]``.

    The function generates a random element of the centralizer
    of ``v_0`` in ``G_x0``.
    """
    a = MM0()
    for i in range(quality):
         pi = rand_pi_mat22()
         x1 = randint(0, 0xfff) & ~0x200
         y1 = randint(0, 0xfff) & ~0x200
         e = randint(0, 2)
         a *= MM0([('p', pi), ('x', x1), ('y', y1), ('l', e)]) 
    return a 



##########################################################################
# Testing function reduce_axis and C function mm_reduce_vector_vp
##########################################################################







def make_axis_testcases():
    #V = V_START.space
    yield V_START.copy(), 1
    yield V("I", 11, 9), 1
    yield V("I", 11, 9), 0
    for i in range(10):
        yield V_START * MM0("r", "G_x0"), i & 1
    for ax in get_sample_axes().values():
        v0 = ax.v15
        for i in range(20):
            yield v0 * MM0("r", "G_x0"), i & 1
    for quality in range(2,11):
        for i in range(3):
              yield  V_START *  MM0("r", quality), i & 1      


ERR_REDUCE_V_AXIS = "Error in function reduce_v_axis_C, status = %d"


ZERO = np.zeros(1, dtype = np.uint32)

def reduce_v_axis_C(v, std_axis):
    va = mm_vector(15, 2)
    vc, work = va[0].data, va[1].data
    mm_op_copy(15, v.data, vc)    
    r = np.zeros(200, dtype = np.uint32)
    res = mm_reduce_vector_vp(ZERO, vc, std_axis, r, work)
    if not 0 < res <= 200:
        raise ValueError(ERR_REDUCE_V_AXIS, res) 
    g = MM0('a', r[:res])
    axis_found = r[res - 1]
    assert axis_found & 0xfe000000 == 0x84000000
    axis_found &= 0x1ffffff
    return g, axis_found


@pytest.mark.axes
def test_reduce_axis(verbose = 0):
    import_all()
    for i, (v, std_axis) in enumerate(make_axis_testcases()):
        if verbose:
            print("\nTest case", i, ", std_axis =", bool(std_axis))
        r, axis = reduce_axis(v.copy(), std_axis, verbose)
        g = MM0('a', r)
        if std_axis:
            assert v * g == V_START
            assert axis == 0x200
        else:
            v1 = v * g
            assert mm_reduce_2A_axis_type(v1.data) >> 24 == 0x21
            a_g = np.zeros(20, dtype = np.uint32)
            r1 = gen_leech2_reduce_type2(axis, a_g)
            g2 = MM0('a', a_g[:r1])
            v2 = v1 * g2
            assert v2 in  [V_START, V_OPP]
            axis2 = gen_leech2_op_word(axis,  a_g, r1)
            assert axis2 & 0xffffff == 0x200
            assert v_leech2_adjust_sign(v2.data, axis2) == axis2

        g1, axis_found = reduce_v_axis_C(v, std_axis)
        assert g1 == g
        assert axis == axis_found   




##########################################################################
# Testing function reduce_baby_axis and mm_reduce_v_baby_axis
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
         a *= rand_Co2()
         a *= MM0('t', 'r')
    return a 

ERR_REDUCE_MM_V_BABY_AXIS = ("Error in function mm_reduce_v_baby_axis"
    ", status = %d")

def mm_reduce_v_baby_axis_C(v, axis = V_START, mode = 0):
    va = mm_vector(15, 2)
    vc, work = va[0].data, va[1].data
    mm_op_copy(15, v.data, vc)    
    r = np.zeros(200, dtype = np.uint32)
    res = mm_reduce_vector_shortcut(1, mode, axis, r)
    if res < 0:
        raise ValueError(ERR_REDUCE_MM_V_BABY_AXIS % res)
    res = mm_reduce_vector_vm(ZERO, vc, r, work)
    if res < 0:
        raise ValueError(ERR_REDUCE_MM_V_BABY_AXIS % res)
    g = MM0('a', r[:res])
    axis_found = r[res - 1]
    assert axis_found & 0xfe000000 == 0x86000000
    axis_found &= 0x1ffffff
    return g, axis_found




def make_baby_testcases():
    def make_pair(v, randomize):
        v, axis = v.copy(), v_start
        if not randomize:
            return v, axis
        g = MM0('r', 'G_x0')
        return v * g, (XLeech2(v_start) * g).ord
               

    yield make_pair(V_OPP, False)
    for i in range(10):
        yield make_pair(V_OPP.copy() * rand_Co2(), i & 1)
    for ax in get_baby_sample_axes().values():
        v0 = ax.v15
        for i in range(5):
            yield make_pair(v0 * rand_BM(), i & 1)
    for quality in range(2,11):
        for i in range(5):
              yield  make_pair(V_OPP.copy() * rand_BM(quality), i & 1)      



@pytest.mark.axes
def test_reduce_baby_axis(verbose = 0):
    import_all()
    for i, (v, axis) in enumerate(make_baby_testcases()):
        if verbose:
            print("\nTest case", i, ", axis =",  hex(axis))
        r, img_axis = reduce_baby_axis(v.copy(), axis, verbose)
        g = MM0('a', r)
        if verbose:
            print("g =",  g)
        v_g = v * g
        if axis == v_start:
            assert V_START * g == V_START
            assert axis == img_axis
            assert v_g == V_OPP
        assert XLeech2(axis) * g == XLeech2(img_axis)
        vt = mm_reduce_2A_axis_type(v_g.data)
        vt = v_leech2_adjust_sign(v_g.data, vt)
        assert  img_axis ^ vt == 0x1000000, (hex(img_axis), hex(vt))


        #vr1 = np.zeros(200, dtype = np.uint32)
        #len_r1 = mm_reduce_v_baby_axis(v.copy().data, axis, vr1)
        #assert len_r1 >= 0, len_r1
        #g1 = MM0('a', vr1[:len_r1])
        #g1 = mm_reduce_v_baby_axis(v.data, axis)
        g1, axis_found =  mm_reduce_v_baby_axis_C(v, axis)

        ok =  g1 == g 

        #print(type(g1), type(g))
        if verbose or not ok:
             vt = mm_reduce_2A_axis_type(v.data) >> 24
             vA = eval_A_vstart(v.data)
             print("Type(v) = 0x%x, value(A) = %d" % (vt, vA))
             print("Op:  ", [hex(x) for x in g.mmdata] )
             print("Fast:", [hex(x) for x in g1.mmdata] ) 
             if not ok:
                 err = "Function mm_reduce_v_baby_axis failed"
                 raise ValueError(err)


##########################################################################
# Testing C function mm_reduce_G_x0
##########################################################################


reduce_time = None

ERR_GX0 = "Error in function mm_reduce_G_x0, status = %d"

def reduce_G_x0(g, mode = 0):
    global reduce_time
    va = mm_vector(15, 2)
    v, work = va[0].data, va[1].data
    r = np.zeros(256, dtype = np.uint32)
    g1, n = g.mmdata, len(g.mmdata)

    t_start = time.perf_counter()
    mm_reduce_load_axis(v, 0);
    res = mm_op_word(15, v, g1, n, 1, work)
    if res < 0: raise ValueError(ERR_GX0 % res)
    res = mm_reduce_vector_vp(ZERO, v, mode, r, work)    
    if res < 0: raise ValueError(ERR_GX0 % res)

    mm_reduce_load_axis(v, 1);
    res = mm_op_word(15, v, g1, n, 1, work)
    if res < 0: raise ValueError(ERR_GX0 % res)
    res = mm_reduce_vector_vm(ZERO, v, r, work)  
    if res < 0: raise ValueError(ERR_GX0 % res)
    reduce_time = time.perf_counter() - t_start
    return MM0('a', r[:res])




def g_complexity(g):
    s = [x for x in g.mmdata if x & 0x70000000 == 0x50000000]
    return len(s)


def make_reduce_testcases():
    for quality in range(1,17):
        for i in range(3):
              yield  MM0('r', quality) , i & 1     


@pytest.mark.axes
def test_reduce_G_x0(verbose = 0):
    import_all()
    for i, (g, mode) in enumerate(make_reduce_testcases()):
        if verbose:
            print("\nTest case", i, "mode =", mode)
        g1 = reduce_G_x0(g, mode)
        g2 = g * g1
        ok = g2.in_G_x0()

        if verbose or not ok:
            print("g =", g)
            print("g1 =", g1 ) 
            print("g2 =", g2 ) 
            print("Run time: %.2f ms, Complexity: %d" % 
                (1000 * reduce_time, g_complexity(g))) 
            if not ok:
                err = "Function mm_reduce_G_x0 failed"
                raise ValueError(err)



