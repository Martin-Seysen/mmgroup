r"""Investigate certain orbits of 2A axes of the monster group

The functions in this module are used to find the orbits of 2A axes
of the subgroup :math:`G_{x0}` (of structure 
:math:`2^{1+24}.\mbox{Co}_1`) of the monster group. According to
|Nor98| there are 12 such orbits.

"""

import sys
import os
import time
import operator
from numbers import Integral
from copy import deepcopy
from math import floor, ceil
from random import randint, shuffle, sample
from collections import defaultdict, OrderedDict
from collections.abc import Iterable
import numpy as np
from operator import __or__
from functools import reduce

if __name__ == "__main__":
    sys.path.append("../../../")

from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2, Parity, PLoop
from mmgroup import mat24, mmv_scalprod
from mmgroup.mat24 import perm_from_dodecads
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup import mat24, GcVector, AutPL, Parity, GCode
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech3to2
from mmgroup.generators import gen_leech3to2_short
from mmgroup.generators import gen_leech2to3_abs
from mmgroup.generators import gen_leech3_add
from mmgroup.generators import gen_leech3_neg
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.generators import gen_leech3_reduce_leech_mod3
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.clifford12 import bitmatrix64_solve_equation
from mmgroup.clifford12 import leech3matrix_kernel_image
from mmgroup.clifford12 import leech3matrix_compress
from mmgroup.clifford12 import leech3matrix_sub_diag
from mmgroup.clifford12 import leech3matrix_vmul
from mmgroup.clifford12 import leech2_matrix_basis
from mmgroup.clifford12 import leech2_matrix_radical
from mmgroup.clifford12 import leech2_matrix_expand
from mmgroup.clifford12 import bitmatrix64_solve_equation
from mmgroup.mm_op import mm_op_load_leech3matrix
from mmgroup.mm_op import mm_op_eval_X_find_abs
from mmgroup.mm_op import mm_op_norm_A
from mmgroup.mm_op import mm_op_eval_A
from mmgroup.mm_op import mm_op_t_A
from mmgroup.mm_op import mm_op_eval_A_rank_mod3

from mmgroup.tests.axes.axis import G, Axis
from mmgroup.tests.axes.axis import G_CENTRAL, G_AXIS, G_AXIS_OPP
from mmgroup.tests.axes.axis import V_AXIS, V_AXIS_OPP

V15 = MMV(15)

#################################################################
# Auxiliary functions for beautifying an axis in Co_1
#################################################################



def find_guide(class_, axis):
    if class_ == '4B':
        Set_1 = axis.find_short(1, 2)
        v0 = Set_1[0]
        for v in Set_1:
            if gen_leech2_type(v ^ v0) == 4:
                return v ^ v0
    if class_ == '12C':
        S1 = set(axis.find_short(1, 2))
        S7 = set(axis.find_short(7, 2))
        short = (S1 & S7) - set([0])
        assert len(short) == 1
        assert len(S7) == 256
        v0 = list(short)[0]
        for v in S7:
            if gen_leech2_type(v) == 2 and gen_leech2_type(v ^ v0) == 4:
                return v ^ v0
    if class_ == "10B":
        #a = np.zeros(3 * 24, dtype = np.uint64)
        #mm_op_load_leech3matrix(axis.p, axis.data, a)
        a = axis.leech3matrix()
        leech3matrix_compress(a, a)
        rad = axis.find_short(4, 2)
        for v in rad:
            if gen_leech2_type(v) == 2:
                v3 = gen_leech2to3_abs(v)
                if v3 > 0:
                   w3 = leech3matrix_vmul(v3, a)
                   #print("w3 =", hex(w3))
                   if w3 > 0:  
                      w = gen_leech3to2(w3)
                      if w > 0 and gen_leech2_type(w) == 4:
                          return w
    if class_ == "6C":
        rad = axis.find_short(3, 2)
        #print([gen_leech2_type(x) for x in rad])
        for v in rad:
            if gen_leech2_type(v) == 4:
                #print("v=", v)
                return v
    return 0



def reduce_leech_mod_3(v):
    g = np.zeros(12, dtype = np.uint32)
    len_g = gen_leech3_reduce_leech_mod3(v, g) >> 48
    assert len_g >= 0
    v_reduced = gen_leech3_op_vector_word(v, g, len_g)
    # print("v_reduced", hex(v_reduced))
    return G('a', g[:len_g])

def get_v3_case_4A(axis):
    ker, isect, _, a_ker = axis.kernel_A(0)
    assert ker == 1
    return leech3matrix_vmul(1, a_ker[23:]) 

def get_v3_case_6A(axis):
    rad  = axis.find_short(5, radical = 1)
    a = np.array(list(map(gen_leech2to3_abs, rad)), dtype = np.uint64)
    assert len(a) == 2
    for v in [gen_leech3_add(a[0], a[1]),
                 gen_leech3_add(a[0], gen_leech3_neg(a[1]))]:
        type_v = gen_leech3to2(v) >> 24
        if not 2 <= type_v <= 4:
            return v
    raise ValueError("Reduction for axis orbit 6A failed")

def get_v3_case_10A(axis):
    ker, isect, _, a_ker = axis.kernel_A(0)
    assert ker == 2
    a = a_ker[22:]
    vl = [leech3matrix_vmul(v0, a) for v0 in [1, 2, 3, 0x1000002]]
    for v in vl:
        typev, v2 =  divmod(gen_leech3to2(v), 1 << 24)
        if not 2 <= typev <= 4:
           return v
    raise ValueError("Reduction for axis orbit 10A failed")

ORBIT_MOD3_CASES = {
    '4A' : get_v3_case_4A,
    '6A' : get_v3_case_6A,
    '10A' : get_v3_case_10A,
}

def do_get_v3_case_2A(axis):
    if len(axis.g.mmdata) == 0:
        return axis
    ker, isect, a_img, _ = axis.kernel_A(0)
    assert ker == 23
    v3 = leech3matrix_vmul(1, a_img)
    typev, v2 =  divmod(gen_leech3to2(v3), 1 << 24)
    assert typev == 2
    g = np.zeros(6, dtype = np.uint32)
    len_g = gen_leech2_reduce_type2(v2, g)
    assert 0 <= len_g <= 6
    axis *= G('a', g[:len_g])
    if (int(axis['B', 2, 3]) + 2) % 15:
        axis *= G('x', 0x200)
    assert axis.v15 ==  axis.v_axis15 
    return axis

#################################################################
# Final permutations and sign changes for beautifying an axis
#################################################################


def solve_gcode_diag(l, generator = 'y'):
    """Solve cocode equation

    Here ``l`` is a list of tupeles ``(i0, i1, k)``.  For an
    unknown Golay code word ``x``, each tuple means an equation
    ``<x, Cocode([i0,i1])> = k``, where ``<.,.>`` is the scalar
    product. If a solution ``x`` exists then the function
    returns the element ``MM0('y', x)`` of the group ``AutPL``
    as an instance of class MM0.
    """ 
    a = np.zeros(len(l), dtype = np.uint64)
    for i, (i0, i1, k) in enumerate(l):
        a[i] = Cocode([i0,i1]).ord + ((int(k) & 1) << 12)
    v = bitmatrix64_solve_equation(a, len(l), 12)
    if v < 0:
        err = "Off-diagonal matrix equation has no solution"
        raise ValueError(err)    
    result = PLoop(v)
    for i0, i1, k in l:
        c = Cocode([i0,i1])
        check = hex(result.ord), hex(c.ord), k
        assert result & c == Parity(int(k)), check
    if generator == 'xy':
        return G('x', result) * G('y', result)
    return G(generator, result)



def permutation(image, preimage = None):
    """Compute permutation mapping ``preimage`` to ``image``

    The function computes an element of the Mathieu group
    :math:`M_{24}` that maps the ``preimage`` to the ``image``.
    That element of the Mathieu group need not be unique
    These two parameters must be list of integers. In case
    ``preimage == None```we take the preimage ``range(24)``.
  
    If such a permutation ``pi`` exists then the function returns
    the premiage ``MM0('p', AutPL(0, pi))`` of ``pi`` in the
    group ``AutPL`` as an instance of class MM0.  
    """
    if preimage is None:
        return G('p', AutPL(0, image))
    else:
        return G('p', AutPL(0, zip(preimage, image), unique=False))

SWAP = G('d', 0x800)

def cond_swap_BC(condition):
    return  SWAP if condition else G()


L12C = [8, 9, 16, 18, 21, 22]  # [2,3] + L12C  is an octad
DODECAD_12C = [0, 1, 4, 5, 6, 7] + L12C
for s in (L12C, DODECAD_12C[:6]):
    assert mat24.vect_to_octad(sum(1 << i for i in [2,3] + s))


OCTAD_10B = [2, 3, 0, 1, 4, 5]

NEG_DODECAD_8F = [12, 16, 20, 1, 2, 3, 5, 6, 7, 9, 10, 11]
DODECAD_8F = [0, 4, 8, 13, 14, 15, 17, 18, 19, 21, 22, 23]


def debug_show_swapped(axis, swap = False, verbose = 1, submatrix = None):
   if submatrix is None:
       from_, to_ = 0, 24
   else: 
       from_, to_ = submatrix
   if swap == 'r':
      swap = randint(0, 1)
   if bool(swap):
       print("Swapping parts B and C for debugging!")
       axis *= cond_swap_BC(bool(swap))
   if verbose:
       print("type axis * t: %s, axis * t^2: %s" % 
             (axis.axis_type(1) , axis.axis_type(2)) )
       for e in range (3):
            print("Orbit of axis * t**%d is %s" % (e, axis.axis_type(e)))
            display_A((axis * G('t', e))['A'])


REF_CLASS = None


def find_perm_6F(i):
    syn = mat24.lsbit24(mat24.syndrome(0x10e | (1 << i)))
    cosyn = mat24.lsbit24(0xf1 & ~(1 << syn) )
    #print("SYN", i, syn, cosyn)
    src = [syn, 1, 2, 3, 8, i, cosyn]
    #print(src)
    return  [0, 1, 2, 3, 8, 9, 4], src 


def postpermute(class_, axis):
    if class_ == "12C":
        # display_A(axis['A'] ) 
        diag = np.diagonal(axis['A'])  % 3
        # print(diag)
        preimg = [i for i, x in enumerate(diag[:8]) if x == 2]
        preimg += [int(i) + 8 for i, x in enumerate(diag[8:]) if x == 1]
        # print("preimage bits", GCode(preimg).bit_list)
        # print("DODECAD", GCode(DODECAD_12C).bit_list)
        pi = perm_from_dodecads(preimg, DODECAD_12C)
        axis *= permutation(pi)
        a = axis['A'] % 3
        l = [(2, 3, a[2,3] != 2)]
        L12C = [8, 9, 16, 18, 21, 22] # [2,3] + L12C  is an octad
        for j in L12C[1:]:
            l += [(L12C[0], j, a[L12C[0], j] != 1)]
        axis *= solve_gcode_diag(l, 'y')
        b = axis['B'] % 3
        l = []
        for j in [0,1,4,5,6,7]:
            l += [(2, j, b[2,j] != 2)]
        for j in L12C[1:]:
            l += [(2, j, b[2,j] != 2)]
        axis *= solve_gcode_diag(l, 'x')
        axis *= cond_swap_BC(axis['B', 0, 8] != axis['C', 0, 8])
    if class_ == "10B":
        diag = np.diagonal(axis['A'])  % 3
        o = [None] * 2
        for i, x in enumerate(diag[:8]):
            if x in [1,2]:
                o[x - 1] = i
            else:
                o.append(i)
        axis *= permutation(OCTAD_10B, o[:6])
        #debug_show_swapped(axis, 0, 1)
        c = axis['C'] % 3
        swaptest =  int(c[12,13] * c[12,14] * c[13,14] % 3)
        # print("swaptest", swaptest)
        axis *= cond_swap_BC(swaptest  == 2 )      

        a = axis['A'] % 3
        l = [ (2, 3, (int(a[2,3]) - 2) & 1)]
        l += [(0, j, (int(a[0,j]) - 1) & 1) for j in [1, 4, 5, 6, 7]]
        axis *= solve_gcode_diag(l)
        #display_A(axis['A'] )
        a = axis['A'] % 3
        l = [0, 1, 2, 3]
        l += [j for j in [4, 5, 6, 7] if a[1, j] == 1]
        axis *= permutation(range(6), l)
        #display_A(axis['A'] )
        b = axis['B'] % 3
        l = [(0, j, (int(b[0,j]) - 1) & 1) for j in range(8, 24)]
        axis *= solve_gcode_diag(l, 'x')
        axis *= cond_swap_BC(axis['B', 0, 8] == axis['C', 0, 8])
    if class_ == "10A":
        axis *= cond_swap_BC(axis['B', 0, 9] == axis['C', 0, 9])
        b = axis['B'] % 3
        l = [(0, j, b[0,j] != 1) for j in range(9,24) if  b[0,j]]
        axis *= solve_gcode_diag(l, 'x')
    if class_ == "6F":
        b = axis['B', 8] % 3
        c = axis['C', 8] % 3
        #print(b)
        for i in range(9, 24):
            if b[i]:
                axis *= cond_swap_BC(b[i] == c[i])
                #print(i)
                break
        dest, src = find_perm_6F(i)             
        axis *= permutation(dest, src)
        l = [(i, i+1, axis['C', i, i+1] % 3) for i in range(8, 24, 2)]
        axis *= solve_gcode_diag(l, 'x')
        #axis.display_sym('B', ind = range(8,12), text = 'orbit 6F, axis B')
    if class_ == "6C":
        diag = np.diagonal(axis['A'])  % 3
        head = []; tail = []
        for i, x in enumerate(diag[:8]):
            if x == 0:
                head.append(i)
            else:
                tail.append(i)
        preimage = head + tail[:3] + [8]
        #print(preimage)
        image = [1, 2, 3, 0, 4, 5, 8]
        axis *=permutation(image, preimage)
        #display_A(axis['A'] )
        a = axis['A'] % 3
        l = [(1, i, (int(a[1,i]) - 2) & 1) for i in (2, 3)]
        axis *= solve_gcode_diag(l)
        axis *= cond_swap_BC(axis['B', 0, 8] == axis['C', 0, 8])
        b = axis['B', 0]
        x = sum((1 << i for i in range(8, 24) if b[i] != 9))
        try:
            axis *= G('x', mat24.vect_to_gcode(x))
        except:
            pass
    if class_ == "8B":
        diag = np.diagonal(axis['A'])
        if diag[0] == 9:
            pi = perm_from_dodecads(DODECAD_8F, NEG_DODECAD_8F)
            axis *= permutation(pi)
            diag = np.diagonal(axis['A'])
        position =  list(diag).index(3)
        dodecad = DODECAD_8F[:]
        ind = dodecad.index(position)
        dodecad[0], dodecad[ind] = dodecad[ind], dodecad[0]
        axis *= permutation(perm_from_dodecads(dodecad, DODECAD_8F))
        axis *= cond_swap_BC(axis['B', 0, 2] == axis['C', 0, 2])
        #display_A(axis['B'] )
        b = axis['B', 0]
        l =  [(0, j, b[j] != 1) for j in range(1, 24)
             if b[j] and j != 12]
        axis *= solve_gcode_diag(l, 'x')
    if class_ == "6A":
        diag = np.diagonal(axis['A'])
        position =  list(diag).index(6)
        if position == 2:
            image = [0, 3, 1, 2, 4] 
        elif position == 3:
            image = [0, 2, 3, 1, 4] 
        else:
            image = [0, 1, 2, 3, 4]
        axis *= permutation(image, [0, 1, 2, 3, 4])
        axis *= cond_swap_BC(axis['B', 1, 0] == axis['C', 1, 0])
        b = axis['B', 1] % 3
        l =  [(1, j, b[j] != 1) for j in range(4, 24)]
        axis *= solve_gcode_diag(l, 'x')
    if class_ == "4B":
        diag = np.diagonal(axis['A']) % 3
        head = []; tail = []
        for i, x in enumerate(diag[:8]):
            if x == 1:
                head.append(i)
            else:
                tail.append(i)
        preimage = head + tail[:4]
        image = [2, 3, 0, 1, 4, 5]
        # print(preimage, image)
        axis *= permutation(image, preimage)
        a = axis['A'] % 3
        l = [(2, 3, a[2, 3] == 1)]
        axis *= solve_gcode_diag(l)
        axis *= cond_swap_BC(axis['B', 1, 0] == axis['C', 1, 0])
        l = [(0, j, 0) for j in range(1, 8)]
        b = axis['B', 8] % 3
        l += [(8, j, b[j] != 1) for j in range(9, 24)]
        axis *= solve_gcode_diag(l, 'x')
    if class_ == "4C":
        axis *= cond_swap_BC(axis['B', 0, 8] == axis['C', 0, 8])
        l = [(j, j+1, axis['B', j, j+1] % 3 != 1) for j in range(8, 24, 2)]
        b = axis['B', 8] % 3
        l +=  [(8, j, b[j] == 2) for j in range(8)]
        b = axis['B', 9] % 3
        l +=  [(9, j, b[j] == 2) for j in range(8)]
        b = axis['B', 10] % 3
        l +=  [(10, j, b[j] == 2) for j in range(4)]
        b = axis['B', 12] % 3
        l +=  [(12, j, b[j] == 2) for j in range(2)]
        b = axis['B', 16] % 3
        l +=  [(16, 0, b[0] == 2)]
        axis *= solve_gcode_diag(l, 'x')
    if class_ == "4A":
        axis *= cond_swap_BC(axis['B', 1, 0] == axis['C', 1, 0])
        b = axis['B', 0] // 3
        l = [(0, j, b[j] == 1) for j in range(1, 24)]
        axis *= solve_gcode_diag(l, 'x')
    return axis

    

#################################################################
# Beautifying an axis
#################################################################



USE_CENTRAL_INVOLUTION = [
    '2B', '4B', '4C', '6C', '8B', '6F', '10B', '12C'
]

def beautify_using_central_involution(class_, axis):
    guide = find_guide(class_, axis)
    zi = axis.central_involution() * axis.g_central
    _, h = zi.conjugate_involution_G_x0(guide)
    axis *= h
    return axis



def central_involution_is_ok(class_, axis, verbose = 0):
    g2 = axis.central_involution() * axis.g_central
    cls_, h = g2.conjugate_involution_G_x0()
    g_test = g2**-1 * g2**h
    ok = g_test == G()
    if verbose and not ok:
        print("BAD case:", class_, cls_, g_test)
    return ok

def beautify_axis(axis, rand = 0, check = False, verbose = 0):
    """Try to beautify a 2A axis

    The function multiples an ``axis`` of type ``Axis`` with an
    element ``g`` of the subgroup :math:`G_{x0}` so that a *nice*
    representative of the :math:`G_{x0}` orbit of the axis is obtained.
    There is no canonical *nicest* representative of an orbit. As a
    general idea, we try to minimize the nonzero off-diagonal elements
    of parts 'A', 'B', and 'C' of the axis. Also, we try to put the
    signs of these of-diagonal elements into an orderly shape.

    Let ``z`` be the central involution of :math:`G_{x0}`, ``t`` be
    the involution corresponding to the axis, and ``i`` be the central
    involution in the dihedral group generated by ``t`` and ``z``. Then
    ``z * i`` is an involution in :math:`G_{x0}`; and we also try to
    change the axis so that ``z * i`` is equal to the representative
    of its involution class in :math:`G_{x0}`. These representatives
    are as in method ``conjugate_involution_G_x0`` of
    class ``Xsp2_Co1``.

    We make no attempt to beautify any of the tags 'TXZY' of the axis. 

    The axis is rebased (using method ``rebase`` of class ``Axis``)
    before performing any operation on it. Thus attribute ``g1`` of
    the axis will contain the group operation performed on the axis.

    If argument ``rand`` is True then we multiply the Axis with a
    random element of :math:`G_{x0}` before beautifying it. If argument
    ``check`` is True then we make a (rather soft) check that the
    properties of the involution ``i`` stated above are satisfied.
    """
    assert isinstance(axis, Axis)
    if REF_CLASS and class_ != REF_CLASS:
        # for debugging only
        return None
    class_ = axis.axis_type()
    axis.rebase()
    if rand:
        axis *= G('r', 'G_x0')
    if verbose > 1:
        print("Input:")
        display_A(axis['A'])
    if class_ == "2A":
        axis = do_get_v3_case_2A(axis)
    if class_ in USE_CENTRAL_INVOLUTION:
        beautify_using_central_involution(class_, axis)
        axis = postpermute(class_, axis)
    if class_ in ORBIT_MOD3_CASES:
        v3 = ORBIT_MOD3_CASES[class_](axis)
        axis *= reduce_leech_mod_3(v3) 
        axis = postpermute(class_, axis)
    if check:
        if not central_involution_is_ok(class_, axis):
            if class_ in USE_CENTRAL_INVOLUTION:
                beautify_using_central_involution(class_, axis)
                axis = postpermute(class_, axis)
        central_involution_is_ok(class_, axis, verbose)

    if verbose:
        for e in range (3):
            print("Orbit of axis * t**%s is %s" % (e, axis.axis_type(e)))
            axis.display_sym(e)
        if not rand and mmv_scalprod(axis.v15, axis.v_axis15) == 0:
            print("Axis may also be used for the Baby monster\n")
    return axis


#################################################################
# Main program
#################################################################


if __name__ == "__main__":
    from mmgroup.tests.axes.sample_axes import g_strings, g_classes
    for g in  g_strings:
        beautify_axis(Axis(g), verbose = 1) 




