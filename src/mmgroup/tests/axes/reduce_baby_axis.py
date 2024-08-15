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
from mmgroup.mm_reduce import mm_reduce_find_type4_axis

from mmgroup.tests.axes.axis import Axis, BabyAxis
from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
from mmgroup.tests.axes.equations import solve_Qx0_equations
from mmgroup.tests.axes.beautify_axes import get_v3_case_10A
from mmgroup.tests.axes.beautify_axes import get_v3_case_6A
from mmgroup.tests.axes.beautify_axes import get_v3_case_4A
from mmgroup.tests.axes.beautify_axes import solve_gcode_diag
G = Axis.group


Q_AXIS_OPP = 0x1000200


BABY_AXES = get_baby_sample_axes()



V_PLUS = 0x200
XL2_VPLUS = XLeech2(V_PLUS)




#################################################################
# Auxiliary functions for beautifying an axis in Co_1
#################################################################



def reduce_leech_mod_3(v):
    g = np.zeros(12, dtype = np.uint32)
    len_g = gen_leech3_reduce_leech_mod3(v, g) >> 48
    assert len_g >= 0
    v_reduced = gen_leech3_op_vector_word(v, g, len_g)
    # print("v_reduced", hex(v_reduced))
    return G('a', g[:len_g])


#################################################################
# Permutation
#################################################################

V_1_2 = mat24.vect_to_cocode(6)
V_1_3 = mat24.vect_to_cocode(10)

 
def permutation(image, preimage):
    """Compute permutation mapping ``preimage`` to ``image``

    The function computes an element of the Mathieu group
    :math:`M_{24}` that maps the ``preimage`` to the ``image``.
    That element of the Mathieu group need not be unique
    These two parameters must be list of integers.
  
    If such a permutation ``pi`` exists then the function returns
    the premiage ``MM0('p', AutPL(0, pi))`` of ``pi`` in the
    group ``AutPL`` as an instance of class MM0.  
    """
    im = [2, 3] + list(image) 
    preim = [2, 3] + list(preimage)
    pi = AutPL(0, zip(preim, im), unique=False)
    pi = G('p', pi)
    if XL2_VPLUS * pi != XL2_VPLUS:
        raise(ValueError, "Permutation does not fix axis v^+")
    return pi

def mul_xleech2_g(v, g):
    """Multiply vector in Leech lattice mod 2 with group element

    The function multiplies the (signed) vector ``v`` in the the
    Leech lattice mod 2 with an element ``g`` of the subgroup
    :math:`G_{x0}` of the Monster. Vector ``v`` and the result
    are encoded in *Leech lattice encoding*.
    """
    return gen_leech2_op_word(v, g.mmdata, len(g.mmdata))


def reduce_23(axis, g, pi = None):
    img_v_plus = mul_xleech2_g(V_PLUS, g)
    assert img_v_plus & 0x7ff800 == 0, hex(img_v_plus)
    syn = mat24.cocode_syndrome(img_v_plus)
    b2 = mat24.lsbit24(syn)
    b3 = mat24.lsbit24(syn & ~(1 << b2))
    #print("gggg", g)
    #(Axis(axis) * g).display_sym('A', text = "Non-baby")
    if pi:
        dest = list(pi[1])
        src = [x for x in pi[0] if x not in (b2, b3)][:len(dest)]
    else:
        src = dest = []
    pi1 = AutPL(0, zip([b2, b3] + src, [2, 3] + dest), 0)
    g0 = G('p', pi1)
    img_v_plus = mul_xleech2_g(img_v_plus, g0) 
    g *=  G('p', pi1)
    if img_v_plus & 0x800000:
        g0 = G('y', 0x200)
        img_v_plus = mul_xleech2_g(img_v_plus, g0)
        g *= g0
    if img_v_plus & 0x1000000:
        g0 = G('x', 0x200)
        img_v_plus = mul_xleech2_g(img_v_plus, g0)
        g *= g0
    assert img_v_plus == V_PLUS, hex(img_v_plus)
    axis *= g
    return axis



def make_pi_4b1():
    """Return dictionary for dealing with orbit '4B1'

    For ``8 <= i < 24``, ``i`` even, let ``c(i)`` be the Golay cocode
    word corresponding to the duad ``(i, i+1)``. The function returns
    a dictionary that maps ``c(i)`` to a permutation ``pi(i)`` with
    certain properties listed below. Actually, the values of the
    dictionary are (even) elements of subgroup ``AutPL`` of the
    Monster that map to suitable permutations ``pi(i)``.

    Pemutation ``pi(i)`` maps ``(i, i+1)`` to  ``(2,3)`` and preserves
    the standard trio ``{(8*j, 8*j+1,...,8+j+7) | j = 0, 1, 2}``.
    Furthermore, ``pi(i)`` maps the octad of that trio containing
    ``0`` the the octad containing ``8``.
    """
    DEST_PI_4B1 = {2, 3, 0, 1, 4, 5, 8}
    pi_4b1 = {}
    for i in range(8, 24, 2):
        c = mat24.vect_to_cocode(3 << i)
        src = [i ^ j for j in range(6)] + [0]
        pi = AutPL(0, zip(src, DEST_PI_4B1))
        pi_4b1[c] = G('p', pi)
    return pi_4b1

PI_4B1 = make_pi_4b1() 

def postpermute(orbit, axis):
    G = axis.group
    if orbit == "2B1":
        d = [i for i, x in enumerate(np.diagonal(axis['A']))
            if x == 4 and i not in [2,3]]
        pi = permutation([0,1,4], d[:3])
        axis *= pi
    if orbit == "2B0":
        d = [i for i, x in enumerate(np.diagonal(axis['A'])) if x == 4]
        syn = mat24.syndrome(sum(1 << x for x in d[:3]) + 0xc)
        d3 = mat24.lsbit24(syn & sum(1 << x for x in d))
        ind = d.index(d3)
        d[3], d[ind] = d[ind], d[3]
        pi = permutation(range(8,13), d[:5])
        axis *= pi
    if orbit == "4A1":
        v3 = get_v3_case_4A(axis)
        g3 = reduce_leech_mod_3(v3)
        reduce_23(axis, g3, ([0], [0]))
    if orbit == "4B1":
        #axis.display_sym(text = "pre permute", mod = 15)
        i = axis.central_involution()
        cls_,h = i.conjugate_involution_G_x0(guide = V_PLUS)
        reduce_23(axis, h, pi = (range(8), [0,1,4]))
        d_list = [j for j, x in enumerate(np.diag(axis['A']))
            if x == 10]
        #print("d_list",  d_list)
        neg = 0 if axis['A', d_list[0], d_list[1]] == 13 else 0x800000
        v = mat24.vect_to_cocode((1 << d_list[0]) | (1 << d_list[1]))
        g = G('c', v ^ neg ^ 0x200) ** -1
        reduce_23(axis, g, pi = (range(8), [0,1,4]))
    if orbit == "4C1":
        i = axis.central_involution()
        _, h = i.conjugate_involution_G_x0(guide = V_PLUS)
        img_v_plus = mul_xleech2_g(V_PLUS, h) 
        h *= PI_4B1[img_v_plus & 0x7fffff]
        reduce_23(axis, h, pi = (range(8), [0,1,4]))
    if orbit == "6C1":
        i = axis.central_involution()
        _, h = i.conjugate_involution_G_x0(guide = V_PLUS)
        reduce_23(axis, h, pi = (range(8), [0,1,4]))
        d = [i for i, x in enumerate(np.diagonal(axis['A'])[:8])
            if x == 3 and i not in [2,3]][0]
        src = [2, 3, d]
        src += [x for x in range(8) if x not in src]
        pi = AutPL(0, zip(src[:5], [2, 3, 1, 0, 4]), 0)
        axis *= G('p', pi)
        solve = [(1, j, axis['A', 1, j] != 8) for j in (2,3)]
        axis *= solve_gcode_diag(solve + [(2,3,0)]) 
    return axis     


#################################################################
# Dictionary for reducing an image of vector v+
#################################################################



def prepare_dict_6A1():
    V_1_2 = mat24.vect_to_cocode(6)
    V_1_3 = mat24.vect_to_cocode(10)
    Z = 0x1000000
    swap = G('x', 0x200)
    pi_1_2 = G('p', AutPL(0, zip([1,2],[2,3]), 0)) 
    pi_1_3 = G('p', AutPL(0, zip([1,3],[2,3]), 0)) 
    d = { 
        V_PLUS: G(1), 
        V_PLUS ^ Z: swap,
        V_1_2: pi_1_2, 
        V_1_2 ^ Z: pi_1_2 * swap,
        V_1_3: pi_1_3, 
        V_1_3 ^ Z: pi_1_3 * swap,
    }
    return d



DICT_REDUCE = {
    '6A1' : get_v3_case_6A,
   '10A1' : get_v3_case_10A,
}

DICT_PERM = {
    '6A1' : prepare_dict_6A1()
}
DICT_PERM['10A1'] = DICT_PERM['6A1']


#################################################################
# Reduce axis
#################################################################



def reduce_baby_axis_G_x0(axis, random = False):
    assert isinstance(axis, BabyAxis)
    G = Axis.group
    axis = axis.copy().rebase()
    orbit = axis.axis_type()
    if random:
        axis *= G('r', 'B & G_x0')
    #axis.display_sym(text = "raw")
    t = mm_reduce_find_type4_axis(axis.v15.data, Q_AXIS_OPP)
    ax_type, c = divmod(t, 0x1000000)
    if orbit == '2A0':
        c ^= 0x200
        axis *= G('c', c)**-1
        if axis['C', 2, 3] != BABY_AXES['2A0']['C', 2, 3]:
            axis *= G('d', 0x800)
    elif orbit in DICT_REDUCE:
        v3 = DICT_REDUCE[orbit](axis)
        g3 = reduce_leech_mod_3(v3)
        g3d = g3.mmdata
        img_v_plus = gen_leech2_op_word(V_PLUS, g3d, len(g3d))
        g3 *= DICT_PERM[orbit][img_v_plus]
        axis *= g3
        g = solve_Qx0_equations(orbit, axis)
        axis *= g
    elif orbit != '2A1':
        axis *= G('c', c)**-1
        axis = postpermute(orbit, axis)
        g = solve_Qx0_equations(orbit, axis)
        axis *= g
    return axis.g1


#################################################################
# Check result of reduction
#################################################################



def check_reduced(axis, forced = True, display = False):
    orbit = axis.axis_type()
    ok = axis.v15 == BABY_AXES[orbit].v15
    if ok:
         if display:
             print("Baby Axis of type %s is reduced" % orbit)
         return True
    d = ""
    for tag in "ABCTXZY":
        if (axis[tag] != BABY_AXES[orbit][tag]).any():
            d += tag
    s = "Baby Axis of type %s is not reduced. Unequal tags: %s"
    s = s % (orbit, d)
    if forced:
        raise ValueError(s)
    if display:
        print(s)
    return False

#################################################################
# Main program
#################################################################



if __name__ == "__main__":
    NTESTS = 20
    DONE_CASES = [
    '2A0', '2A1', '2B1', '2B0', '10A1', '6A1', '4A1', '4B1', '4C1', '6C1'
    ]
    ACTIVE_CASES = []
    for orbit, axis in BABY_AXES.items():
      print("Orbit", orbit)
      for i in range(NTESTS):
         if orbit not in DONE_CASES + ACTIVE_CASES:
             continue
         axis1 = axis * axis.group('r', 'B & G_x0')
         g = reduce_baby_axis_G_x0(axis1)
         axis2 = axis1 * g
         if orbit in DONE_CASES:
              check_reduced(axis2, forced = True)
              continue
         forced = orbit in DONE_CASES
         ok = check_reduced(axis2, forced = forced, display = True)
         if not ok:
             for tag in  "ABC": # (0,1,2):
                 axis1.display_sym(tag, mod = 15)
                 axis.display_sym(tag, mod = 15, text = "Ref " + tag)
             #axis.display_sym('A', mod = 15, text = "Sample")


