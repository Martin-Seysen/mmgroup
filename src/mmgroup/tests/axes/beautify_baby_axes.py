r"""Investigate certain orbits of 2A axes of the monster group

The functions in this module are used to find the orbits of 2A axes
of the subgroup :math:`G_{x0}` (of structure 
:math:`2^{1+24}.\mbox{Co}_1`) of the monster group. According to
|Nor98| there are 12 such orbits.

"""

import sys
import os
import time
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
from mmgroup import mat24, Xsp2_Co1
from mmgroup.mat24 import perm_from_dodecads
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup import mat24, GcVector, AutPL, Parity, GCode
from mmgroup.generators import gen_leech2_op_atom, gen_leech2_op_word
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_reduce_type2_ortho
from mmgroup.generators import gen_leech3to2
from mmgroup.generators import gen_leech2to3_abs
from mmgroup.generators import gen_leech3_add
from mmgroup.generators import gen_leech3_neg
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.generators import gen_leech3_reduce_leech_mod3
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
from mmgroup.mm_op import mm_op_t_A

from mmgroup.tests.axes.axis import G, V15, Axis, BabyAxis
from mmgroup.tests.axes.axis import G_CENTRAL, G_AXIS, G_AXIS_OPP
from mmgroup.tests.axes.axis import g_central, g_axis, g_axis_opp
from mmgroup.tests.axes.axis import V_AXIS, V_AXIS_OPP
from mmgroup.tests.axes.axis import v_axis, v_axis_opp
from mmgroup.tests.axes.axis import v_axis15, v_axis_opp15


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
    if class_ == "6C":
        rad = axis.find_short(3, 2)
        #print([gen_leech2_type(x) for x in rad])
        for v in rad:
            if gen_leech2_type(v) == 4:
                #print("v=", v)
                return v
    return 0


def kernels(axis):
    d = {}
    for i in range(3):
        #a = np.zeros(3 * 24, dtype = np.uint64)
        #mm_op_load_leech3matrix(axis.p, axis.data, a)
        a = axis.leech3matrix()
        leech3matrix_sub_diag(a, i, 0)
        leech3matrix_sub_diag(a, 2, 24)
        x = leech3matrix_kernel_image(a)
        leech3matrix_compress(a, a)
        img, ker, isect = x & 255, (x >> 8) & 255, x >> 16
        assert ker + img == 24
        d[i] = ker, isect, a[:24], a[24:48]
    return d

def reduce_leech_mod_3(v):
    g = np.zeros(12, dtype = np.uint32)
    len_g = gen_leech3_reduce_leech_mod3(v, g) >> 48
    assert len_g >= 0
    v_reduced = gen_leech3_op_vector_word(v, g, len_g)
    # print("v_reduced", hex(v_reduced))
    return G('a', g[:len_g])


def from_subspace_mod3(a, transformed = True):
    if len(a) == 1:
       return leech3matrix_vmul(1, a)
    if len(a) == 2:
        result = 0
        if transformed:
            vl = [leech3matrix_vmul(v0, a) for v0 in [1, 2, 3, 0x1000002]]
        else:
            vl = [a[0], a[1], gen_leech3_add(a[0], a[1]),
                 gen_leech3_add(a[0], gen_leech3_neg(a[1]))]
        for v in vl:
            typev, v2 =  divmod(gen_leech3to2(v), 1 << 24)
            if not 2 <= typev <= 4:
                result = v
        return result

def longest_in_kernel(kernels):
    a = None
    for dim_ker, dim_isect, a_img, a_ker in kernels.values():
        dim_img = 24 - dim_ker
        if 0 < dim_ker <= 2:
            a = a_ker[dim_img : dim_img + dim_ker]
        elif 0 < dim_isect <= 2:
            a = a_img[:dim_isect]
        if a is None:
            continue
        # print("dim", len(a), dim_ker, dim_isect)
        v = from_subspace_mod3(a)
        if v:
            return v 
    return 0;      





def alternative_v3(class_, axis):
    if class_ == "6A":
        rad  = axis.find_short(5, radical = 1)
        a = np.array(list(map(gen_leech2to3_abs, rad)), dtype = np.uint64)
        b = np.array(list(map(gen_leech3to2, a)), dtype = np.uint64)
        return from_subspace_mod3(a, transformed = False)
    return 0


#################################################################
# Final permutations and sign changes for beautifying an axis
#################################################################

# Yet under construction!!!

def solve_gcode_diag(l, generator = 'y'):
    """Solve cocode equation

    Here ``l`` is a list of tupeles ``(i0, i1, k)``.  For an
    unknown Golay code word ``x``, each tuple means an equation
    ``<x, Cocode([i0,i1])> = k``, where ``<.,.>`` is the scalar
    product. If a solution ``x`` exists then the function
    returns the element ``MM0('y', x)`` of the group ``AutPL``
    as an instance of class G.
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



def postpermute(class_, axis):
    if class_ == "10A":
        axis *= cond_swap_BC(axis['B', 0, 9] != axis['C', 0, 9])
        b = axis['B'] % 3
        l = [(0, j, b[0,j] != 1) for j in range(9,24) if  b[0,j]]
        axis *= solve_gcode_diag(l, 'x')
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
        l = [(1, i, (a[1,i] - 2) & 1) for i in (2, 3)]
        axis *= solve_gcode_diag(l)
        axis *= cond_swap_BC(axis['B', 0, 8] != axis['C', 0, 8])
        b = axis['B', 0]
        x = sum((1 << i for i in range(8, 24) if b[i] != 9))
        try:
            axis *= G('x', mat24.vect_to_gcode(x))
        except:
            pass
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
        axis *= cond_swap_BC(axis['B', 1, 0] != axis['C', 1, 0])
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
        l = [(2, 3, (a[2, 3] - 1) & 1)]
        axis *= solve_gcode_diag(l)
        axis *= cond_swap_BC(axis['B', 1, 0] == axis['C', 1, 0])
        l = [(0, j, 0) for j in range(1, 8)]
        b = axis['B', 8] % 3
        l += [(8, j, b[j] != 1) for j in range(9, 24)]
        axis *= solve_gcode_diag(l, 'x')
    if class_ == "4C":
        axis *= cond_swap_BC(axis['B', 0, 8] != axis['C', 0, 8])
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
# Cases 2A and 2B
#################################################################

def beautify_2A0(axis):
    ker, isect, a_img, _ = axis.kernel_A(0)
    assert ker == 23
    v3 = leech3matrix_vmul(1, a_img)
    typev, v2 =  divmod(gen_leech3to2(v3), 1 << 24)
    assert typev == 2
    g = np.zeros(6, dtype = np.uint32)
    len_g = gen_leech2_reduce_type2_ortho(v2, g)
    assert 0 <= len_g <= 6
    axis *= G('a', g[:len_g])
    # Todo: 
    if (axis['C', 2, 3] - 2) % 15:
         # Todo: Change sign of axis!!!!!
         pass
    return axis


#################################################################
# Beautifying an axis
#################################################################


# Yet under construction!!!

def beautify_baby_axis(class_, g, verbose = 0, rand = 0):
    if REF_CLASS and class_ != REF_CLASS:
        # for debugging only
        return None

    axis = BabyAxis(g)
    if class_ == "2A1":
        return axis
    if verbose:
        print("\nOrbit" + class_)
    axis.rebase()
    if rand:
        axis *= G('r', 'G_x0')
    if verbose > 1:
        print("Input:")
        display_A(axis['A'])

    v3 = alternative_v3(class_, axis)
    if v3 == 0:
        d = kernels(axis)
        v3 = longest_in_kernel(d)

        #print("ker", hex(v3))
        #print(class_, [(i, x[:2]) for i, x in d.items()])
    guide = find_guide(class_, axis)
    zi = axis.central_involution(guide) * axis.g_central
    iclass, g2 = zi.conjugate_involution_G_x0(guide)
    if verbose:
        print("guide =", hex(guide))
        print("iclass = ", iclass)
        s = "" if in_Baby(g2) else " NOT"
        print("Central involution is%s in Baby Monster" % s)
    if verbose > 1:
        for e in range (3):
            print("Orbit of axis * t**%d is %s" % (e, axis.axis_type(e)))
            display_A((axis * G('t', e))['A'])
    return axis


#################################################################
# Main program
#################################################################


if __name__ == "__main__":
    from mmgroup.tests.axes.baby_sample_axes import g_strings, g_classes
    for class_, g in zip(g_classes, g_strings):
        beautify_baby_axis(class_, g, verbose = 1) 




