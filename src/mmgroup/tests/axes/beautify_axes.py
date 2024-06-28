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
from mmgroup import mat24
from mmgroup.mat24 import perm_from_dodecads
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup import mat24, GcVector, AutPL, Parity, GCode
from mmgroup.generators import gen_leech2_type
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

from mmgroup.tests.axes.get_sample_axes import G_CENTRAL
from mmgroup.tests.axes.get_sample_axes import G_AXIS
from mmgroup.tests.axes.get_sample_axes import V_AXIS
from mmgroup.tests.axes.get_sample_axes import G_AXIS_OPP
from mmgroup.tests.axes.get_sample_axes import V_AXIS_OPP
from mmgroup.tests.axes.get_sample_axes import g_central
from mmgroup.tests.axes.get_sample_axes import g_axis
from mmgroup.tests.axes.get_sample_axes import g_axis_opp



V15 = MMV(15)



v_axis = MMVectorCRT(16, V_AXIS)
v_axis_opp = MMVectorCRT(16, V_AXIS_OPP)
v_axis15 = V15(V_AXIS)
v_axis_opp15 = V15(V_AXIS_OPP)

#################################################################
# Display a two-dimensional matrix
#################################################################

def display_A(A):
   """Display a two-dimensional integer matrix A"""
   imax, jmax = A.shape
   fmt = [4 if max(abs(A[:,j])) > 99 else 2 for j in range(jmax)]
   for i in range(imax):
      print(" ", end = "")
      for j in range(jmax):
          l = fmt[j]
          if i == j or A[i,j] != 0:
              print("%*d" % (fmt[j], A[i,j]), end = " ")
          else:
              print("%*s" % (fmt[j], "."), end = " ")
      print("")

def sym_part(v, part):
    if isinstance(part, str) and part in 'ABC':
        return v[part]
    if  isinstance(part, Integral):
        a = MMV(v.p)(0)
        mm_op_t_A(v.p, v.data, part, a.data)
        return a['A']
    ERR =  "Cannot display a part of % of a vector"
    raise TypeError(ERR % type(part))

#################################################################
# Classes for modelling an axis
#################################################################


class Axis:
    g_axis_start = g_axis
    v15_start = v_axis15
    g_central = g_central
    constant = False
    ERR1 = "This axis is constant and may not be changed"
    """Models an axis modulo 15"""
    def __init__(self, g = 1):
        self.g0 = MM0()
        self.g1 = MM0(g)
        self.v15 = v_axis15 * self.g
    @property
    def g(self):
        return self.g0 * self.g1
    @property
    def g_axis(self):
        return self.g_axis_start ** self.g
    @property
    def norm_A_mod15(self):
        return mm_op_norm_A(15, self.v15.data)
    def copy(self):
        axis = deepcopy(self)
        if axis.constant:
            del axis.constant
        return axis
    def __imul__(self, g):
        if self.constant:
            raise TypeError(self.ERR1)
        g = MM0(g)
        self.g1 *= g
        self.v15 *= g
        return self
    def __mul__(self, g):
        return self.copy().__imul__(g)
    def __getitem__(self, index):
        return self.v15.__getitem__(index)
    def rebase(self, reduce = False):
        """Rebase group element in axis

        If `´reduce`` is true then
        we try to shorten the group element ``self.g`` with the
        'magic' *mmgroup* functionality.
        This class will be used for proofs; so we'd better check
        the correctness of that magic. We do not change ``self.g``
        if the shortened group element does not map the start
        axis to the current axis.
        """
        if self.constant:
            raise TypeError(self.ERR1)
        self.g0 *= self.g1
        self.g1 = MM0()
        if reduce:
            try:
                from mmgroup import MM
                g0 = MM0(MM(self.g0).reduce())
                assert self.v15_start * g0 == self.v15
                self.g0 = g0
            except:
                import warnings
                W = "Reducing an axis with mmgroup has failed"
                warnings.warn(W, UserWarning)
        return self
    def central_involution(self, guide=0):
        c = g_axis ** self.g * g_central
        _, h = c.half_order()
        result =  h.conjugate_involution_G_x0(guide)
        return result
    def find_short(self, value, radical = 0, verbose = 0):
        return find_short(self.v15, value, radical, verbose)
    def leech3matrix(self):
        a = np.zeros(3 * 24, dtype = np.uint64)
        mm_op_load_leech3matrix(self.v15.p, self.v15.data, a)
        return a
    def axis_type(self, e = 0):
        return self.v15.axis_type(e) 
    def axis_in_space(self, space, *args):
        all_args = args + (V_AXIS,)
        return space(*all_args) * self.g
    def display_sym(self, part = 'A', diff = 0, ind = None, text = "", end = ""):
        if text:
            print(text)
        a =  sym_part(self.v15, part)
        if diff:
            a = (a - sym_part(diff.v15, part)) % 15 
        if ind is not None:
            a = np.array([[a[i, j] for j in ind] for i in ind])
        display_A(a)
        print(end)




def find_short(v, value, radical = 0, verbose = 0):
    """Return certain array of short vectors in Leech lattice mod 2

    Let ``v`` be a vector in the rep of the moster mod 15.

    The function returns the list of all vectors ``v2`` in the Leech
    lattice mod 2 such that the absolute value of the entry of ``v``
    corresponding a unit vector labelled by ``v2`` is equal to
    ``value``. That list is returned as a numpy array.

    If ``radical`` is 1 then the function returns a basis of the
    radical of the space generated by the vectors described above.
    If ``radical`` is 2 then the function returns the whole radical
    as an array of vectors. 
    """
    MAXLEN = 1000
    short = np.zeros(MAXLEN + 1, dtype = np.uint32)
    l = mm_op_eval_X_find_abs(15, v.data, short, MAXLEN+1, value, 0)
    if  l <= 0 or l > MAXLEN:
        return short[:0]
    short = short[:l]
    basis = np.zeros(24, dtype = np.uint64)
    dim = leech2_matrix_basis(short, l, basis, 24)
    if  dim <= 0 or dim >= 24:
        return short[:0]
    if radical or verbose: 
        basis_rad = np.zeros(24, dtype = np.uint64)
        dim_rad = leech2_matrix_radical(short, l, basis_rad, 24)
    if verbose:
        print("The set S_%d has size %3d and dim %2d, radical dim %2d" %
             (value, l, dim, dim_rad))
    if radical == 1:
        return basis_rad[:dim_rad]
    if radical > 1:
        rad = np.zeros(1 << dim_rad, dtype = np.uint32)
        leech2_matrix_expand(basis_rad, dim_rad, rad)
        return rad
    return short



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
    return MM0('a', g[:len_g])


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
        return MM0('x', result) * MM0('y', result)
    return MM0(generator, result)



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
        return MM0('p', AutPL(0, image))
    else:
        return MM0('p', AutPL(0, zip(preimage, image), unique=False))

SWAP = MM0('d', 0x800)

def cond_swap_BC(condition):
    return  SWAP if condition else MM0()


DODECAD_12C = [0, 1, 4, 5, 6, 7, 10, 11, 17, 19, 21, 22]
# print(GCode(DODECAD_12C).bit_list)

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
            display_A((axis * MM0('t', e))['A'])  


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
        preimg += [i + 8 for i, x in enumerate(diag[8:]) if x == 1]
        # print("preimage bits", GCode(preimg).bit_list)
        # print("DODECAD", GCode(DODECAD_12C).bit_list)
        pi = perm_from_dodecads(preimg, DODECAD_12C)
        axis *= permutation(pi)
        a = axis['A'] % 3
        l = [(2, 3, a[2,3] != 2)]
        for j in [11, 17, 19, 21, 22]:
            l += [(10, j, a[10,j] != 1)]
        axis *= solve_gcode_diag(l, 'y')
        b = axis['B'] % 3
        l = []
        for j in [0,1,4,5,6,7]:
            l += [(2, j, b[2,j] != 2)]
        for j in [10, 11, 17, 19, 21]:
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
        swaptest =  c[12,13] * c[12,14] * c[13,14] % 3
        # print("swaptest", swaptest)
        axis *= cond_swap_BC(swaptest  == 2 )      

        a = axis['A'] % 3
        l = [ (2, 3, (a[2,3] - 2) & 1)]
        l += [(0, j, (a[0,j] - 1) & 1) for j in [1, 4, 5, 6, 7]]
        axis *= solve_gcode_diag(l)
        #display_A(axis['A'] )
        a = axis['A'] % 3
        l = [0, 1, 2, 3]
        l += [j for j in [4, 5, 6, 7] if a[1, j] == 1]
        axis *= permutation(range(6), l)
        #display_A(axis['A'] )
        b = axis['B'] % 3
        l = [(0, j, (b[0,j] - 1) & 1) for j in range(8, 24)]
        axis *= solve_gcode_diag(l, 'x')
        axis *= cond_swap_BC(axis['B', 0, 8] != axis['C', 0, 8])
    if class_ == "10A":
        axis *= cond_swap_BC(axis['B', 0, 9] != axis['C', 0, 9])
        b = axis['B'] % 3
        l = [(0, j, b[0,j] != 1) for j in range(9,24) if  b[0,j]]
        axis *= solve_gcode_diag(l, 'x')
    if class_ == "6F":
        b = axis['B', 8] % 3
        c = axis['C', 8] % 3
        #print(b)
        for i in range(9, 24):
            if b[i]:
                axis *= cond_swap_BC(b[i] != c[i])
                #print(i)
                break
        dest, src = find_perm_6F(i)             
        #axis.display_sym('B', ind = range(8,24), text = 'axis B')
        axis *= permutation(dest, src)
        if axis['B', 8, 9] % 3 != 1:
            axis *= solve_gcode_diag([(8, i, 1) for i in range(9, 24, 2)])
        #print(i, [axis.axis_type(e) for e in range(3)])
        #axis.display_sym('B', ind = range(8,24), text = 'axis B')
        #axis.display_sym('C', ind = range(8,24), text = 'axis C')
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
            axis *= MM0('x', mat24.vect_to_gcode(x))
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
        axis *= cond_swap_BC(axis['B', 0, 2] != axis['C', 0, 2])
        #display_A(axis['B'] )
        b = axis['B', 0]
        l =  [(0, j, b[j] != 1) for j in range(1, 20) if b[j]]
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

    return axis

    

#################################################################
# Beautifying an axis
#################################################################


def beautify_axis(class_, g, verbose = 0, rand = 0):
    if REF_CLASS and class_ != REF_CLASS:
        # for debugging only
        return None

    axis = Axis(g)
    if class_ == "2A":
        return axis
    if verbose:
        print("\nOrbit" + class_)
    axis.rebase()
    if rand:
        axis *= MM0('r', 'G_x0')
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
    #print("guide =", hex(guide))
    iclass, g2 = axis.central_involution(guide)
    if v3 == 0 or 'o' in iclass:
        axis *= g2
        #print(iclass)
    if v3 and class_ not in ["6C", "6F", "12C"]:
        axis *= reduce_leech_mod_3(v3)
    axis = postpermute(class_, axis)
    #print(axis.g)
    if verbose:
        for e in range (3):
            print("Orbit of axis * t**%d is %s" % (e, axis.axis_type(e)))
            display_A((axis * MM0('t', e))['A'])
    return axis


#################################################################
# Main program
#################################################################


if __name__ == "__main__":
    from mmgroup.tests.axes.sample_axes import g_strings, g_classes
    for class_, g in zip(g_classes, g_strings):
        beautify_axis(class_, g, verbose = 1) 




