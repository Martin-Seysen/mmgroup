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
from mmgroup import mat24, mmv_scalprod
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


MM_INITIALIZED = False
G = MM0
G_SHORTEN = False
MM = None

def set_axis_group(group = None, shorten = True):
    """Set group acting on instances of class ``Axis``.

    The ordinary user almost never has a reson to call this
    function. Here parameter ``group`` should be ``None``
    (default), or  one of classes ``MM`` or ``MM0`` defined 
    im module ``mmgroup`. By default we use class ``MM`` if
    present and class ``MM0`` otherwise.

    More details about use cases for class ``MM0`` are yet
    to be documented!
    """ 
    global G, G_SHORTEN, MM, MM_INITIALIZED, Axis
    if shorten or group != MM0:
        try:
            from mmgroup import mm_reduce
            from mmgroup import MM as original_MM
            G_SHORTEN = True
            MM = original_MM
            G = MM0 if group == MM0 else original_MM
        except:
            if group not in [None, MM0]:
                ERR = "Class MM is not implemented"
                raise NotImplementedError(ERR)
            G = MM0
            G_SHORTEN = False
            MM = None
    Axis.group = G
    MM_INITIALIZED = True


#################################################################
# Display a two-dimensional matrix
#################################################################

def display_A(A):
   """Display a two-dimensional integer matrix A"""
   imax, jmax = A.shape
   fmt = []
   for j in range(jmax):
       lengths = [len(str(int(A[i,j]))) for i in range(imax)]
       fmt.append(max(2, max(lengths)))
   #fmt = max( (A[:,j]) for j in range(jmax)])
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
        if isinstance(v, MMVector):
            a = MMV(v.p)(0)
            mm_op_t_A(v.p, v.data, part, a.data)
        else:
            a = deepcopy(v) * G('t', part)
        return a['A']
    ERR =  "Cannot display a part of % of a vector"
    raise TypeError(ERR % type(part))

#################################################################
# Classes for modelling an axis
#################################################################


class Axis:
    """Model a 2A axis in the representation of the Monster mod 15

    An instance of this class is an arbitrary 2A axis ``a`` of the
    Monster group with entries taken mod 15. Attribute ``a.v15`` of
    type ``|MMVector|`` is the vector in the representation of the
    Monster corrsponding to the axis. Apart from that vector we also
    maintain a property ``a.g`` of class ``MM0`` such that the axis
    is ``a.v15`` equal to the product ``a.v15_start * a.g``. Here
    ``a.v15_start`` is the standard 2A axis :math:`v^+` defined in
    :cite:`Seysen22`. One may multiply an axis with an element of
    the Monster group. One may abbreviate ``a.v15[x]`` to ``a[x]``
    for any value ``x``.

    This class is also used for proofs. Therefore compute in the
    Monster using instances of class ``MM0``. That class implements
    operations of the Monster defined in :cite:`Seysen20` only.
    Note that this way we cannot shorten the word ``a.g``
    representing an element of the Monster. On may use method
    ``rebase`` for shortening ``a.g``.

    The constructor may of this class may be:

    * An element :math:`g` of the Monster. Then we constuct the
      axis c:math:`v^+ * g`.

    * The string ``'i'`` followed by a 2A involution :math:`g` in
      the Monster. Then we construct the axis corresponding
      to :math:`g`.

    * An instance of this class. Then we construct a copy of that
      instance.
    """
    v15_start = v_axis15
    v15_start_name = V_AXIS 
    constant = False
    group = G
    ERR1 = "This axis is constant and may not be changed"
    def __init__(self, g = 1, invol = None):
        if isinstance(g, Axis):
            self.g1 = G(g.g)
        elif g == 'i':
            from mmgroup import MM
            g1 = MM(invol)
            i, h = g1.conjugate_involution()
            if i != 1:
                ERR = "Element is not a 2A involution in the Monster"
                raise ValueError(ERR)
            self.g1 = G(h**-1)
        else:
            self.g1 = G(g)
        self.g0 = G()
        self.v15 = v_axis15 * self.g
    @property
    def g_central(self):
        """Return central involution of the group ``G_x0``"""
        return G(G_CENTRAL)
    @property
    def g_axis_start(self):
        """Return the fixed 2A axis ``v^+``"""
        return G(G_AXIS)
    @property
    def g(self):
        """Return ``g`` with ``a.v15 = a.v15_start * a.g``"""
        return self.g0 * self.g1
    @property
    def g_axis(self):
        """Return the 2A involution corresponding to the 2A axis"""
        return self.g_axis_start ** self.g
    @property
    def group(self):
        """Return class implementing the Monster group acting on axes

        This is usually class |MM|, and need not be changed. One may
        change this by calling function ``set_axis_group``.
        """
        return G
    @property
    def norm_A_mod15(self):
        """Return norm of part A of the axis (modulo 15)."""
        return mm_op_norm_A(15, self.v15.data)
    def copy(self):
        """Return a (deep) copy of the axis"""
        axis = deepcopy(self)
        if axis.constant:
            del axis.constant
        return axis
    def __imul__(self, g):
        if self.constant:
            raise TypeError(self.ERR1)
        g = G(g)
        self.g1 *= g
        self.v15 *= g
        return self
    def __mul__(self, g):
        return self.copy().__imul__(g)
    def __getitem__(self, index):
        return self.v15.__getitem__(index)
    def rebase(self):
        """Rebase group element in axis

        For an instance ``a`` of class ``Axis`` we try to shorten the
        group element ``a.g`` with the full *mmgroup* functionality.
        Here we may also change ``a.g`` to a shorter element mapping
        the standard axis to the current axis.

        Class ``Axis`` will be used for proofs; so we always check
        the correctness of such operations with the methods of
        class ``MM0``.

        We store ``a.g`` as a product ``a.g1 * a.g2``. Construction of
        ``a`` and multiplication of ``a`` with a group element affects
        ``a.g2`` only. This method puts ``a.g1 = a.g`` and ``a.g2 = 1``
        and reduces ``a.g1`` using the features of class |MM|.
        """
        if self.constant:
            raise TypeError(self.ERR1)
        self.g0 *= self.g1
        self.g1 = G()
        try:
            from mmgroup import MM
            from mmgroup.mm_reduce import mm_reduce_vector_vp 
            v0 = np.zeros(1, dtype = np.uint32)
            v = self.v15.copy()
            w = V15(0)
            g = np.zeros(128, dtype = np.uint32)
            lg = mm_reduce_vector_vp(v0, v.data, 1, g, w.data)
            assert 0 <= lg < 128
            g0 = G(MM('a', g[:lg]) ** -1)
            assert self.v15_start * g0 == self.v15
            self.g0 = g0
        except:
            import warnings
            W = "Reducing an axis with mmgroup has failed"
            warnings.warn(W, UserWarning)
            raise
        return self
    def central_involution(self, guide=0):
        """Return central involution of dihedral group

        The method returns the central involution of the dihedral
        group generated by ``a.g_axis`` and ``a.g_central``. For an
        instance ``a`` of this class, the group element ``a.g_axis``
        is the 2A involution corresponding to the axis, and
        ``a.g_central`` the central involution in :math:`G_{x0}`.
        """
        c = self.g_axis * self.g_central
        _, h = c.half_order()
        return h
    def find_short(self, value, radical = 0, verbose = 0):
        """Return array of short vectors in Leech lattice mod 2

        Let ``v = a.v15`` be the axis given by instance ``a`` of this
        class.  The function returns the list of all vectors ``v2`` in
        the Leech lattice mod 2 such that the absolute value of the
        entry of ``v`` corresponding a unit vector labelled by ``v2``
        is equal to ``value``. That list is returned as a numpy array.

        If ``radical`` is 1 then the function returns a basis of the
        radical of the space generated by the vectors described above.
        If ``radical`` is 2 then the function returns the whole
        radical as an array of vectors.
        """
        return find_short(self.v15, value, radical, verbose)
    def leech3matrix(self):
        """Return the 'A' part of the axis.

        The function returns the part of with tag 'A' of the axis
        as a matrix. That matrix a will be given in matrix mod 3
        encoding as documented in the header of
        file ``leech3matrix.c``.
        """
        a = np.zeros(3 * 24, dtype = np.uint64)
        mm_op_load_leech3matrix(self.v15.p, self.v15.data, a)
        return a
    def axis_type(self, e = 0):
        """``a.axis_type(e)`` is equivalent to ``a.v15.axis_type(e)``"""
        return self.v15.axis_type(e) 
    def axis_in_space(self, space, *args):
        """Return axis as an element of a ``space``

        The function returns the axis as an element of the space given
        by ``space(*args)``.
        """
        all_args = args + (self.v15_start_name,)
        return space(*all_args) * self.g
    def display_sym(self, part = 'A', diff = 0, ind = None, text = "", end = ""):
        """Display part of axis as a symmetric 24 times 24 matrix

        if ``part`` is 'A', 'B', or 'C' then the matrix corresponding to
        that part (modulo 15) is displayed. If ``part`` is an integer
        then the 'A' part of the matrix (multiplied with group element
        ``MM('t', part)``) is displayed.

        If ``diff`` is is an instance of class ``Axis`` then ``diff``
        is subtracted from the axis. Parameter ``ind`` is an optional
        subset of the indices of entries of part ``A`` of the axis.

        ``text`` is an optional text to be displayed; ``end`` is an
        optional end mark.
        """
        if text:
            print(text)
        a =  sym_part(self.v15, part)
        if diff:
            a = (a - sym_part(diff.v15, part)) % 15 
        if ind is not None:
            a = np.array([[a[i, j] for j in ind] for i in ind])
        display_A(a)
        print(end)
    def display_sym_integral(self, part = 'A', text = "", end = ""):
        """Display part of 256 * axis as a symmetric 24 times 24 matrix

        Operaration and parameters are as in method ``display_sym``,
        but the exact symmetric matrix (multiplied with 256)  is
        displayed instead.
        """
        if text:
            print(text)
        v = self.axis_in_space(MMVectorCRT, 20) 
        a = 256 * sym_part(v, part)
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




if not MM_INITIALIZED:
    set_axis_group()


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



def kernel_A(axis, d):
    """Compute invariants of a matrix related to part 'A' of an axis.

    Let 'A' be the symmetric 24 times 24 matrix corresponding to part
    'A' of an axis, with entries taken modulo 3.

    Put 'M' = 'A' - 'd' * 'E', where  'd' is a integer and 'E' is the
    unit matrix. We return a tuple '(ker, isect, M_img, M_ker)', where
    'ker' is the dimension of the kernel of 'M', 'isect' is the
    dimension of the intersection of the kernel and the image
    of 'M'; and 'M_img' and 'M_ker' are 24 times 24 matrices. Row 'i'
    of matrix 'M_img' is the image of row 'i' of 'M_ker' under the
    symmetric matrix 'M'. The first 'isect' rows of 'M_img' and
    'M_ker' are equal; they are a basis of the intersection of the
    kernel and the image of 'M'. The first '24 - ker' rows of 'M_img'
    are a basis of the image of 'M' The last 'ker' rows of 'M_ker'
    are a basis of the kernel of 'M'.
    """
    a = axis.leech3matrix()
    leech3matrix_sub_diag(a, d, 0)
    leech3matrix_sub_diag(a, 2, 24)
    x = leech3matrix_kernel_image(a)
    leech3matrix_compress(a, a)
    img, ker, isect = x & 255, (x >> 8) & 255, x >> 16
    assert ker + img == 24
    return ker, isect, a[:24], a[24:48]


def reduce_leech_mod_3(v):
    g = np.zeros(12, dtype = np.uint32)
    len_g = gen_leech3_reduce_leech_mod3(v, g) >> 48
    assert len_g >= 0
    v_reduced = gen_leech3_op_vector_word(v, g, len_g)
    # print("v_reduced", hex(v_reduced))
    return G('a', g[:len_g])

def get_v3_case_4A(axis):
    ker, isect, _, a_ker = kernel_A(axis, 0)
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
    ker, isect, _, a_ker = kernel_A(axis, 0)
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
    ker, isect, a_img, _ = kernel_A(axis, 0)
    assert ker == 23
    v3 = leech3matrix_vmul(1, a_img)
    typev, v2 =  divmod(gen_leech3to2(v3), 1 << 24)
    assert typev == 2
    g = np.zeros(6, dtype = np.uint32)
    len_g = gen_leech2_reduce_type2(v2, g)
    assert 0 <= len_g <= 6
    axis *= G('a', g[:len_g])
    if (axis['B', 2, 3] + 2) % 15:
        axis *= G('x', 0x200)
    assert axis.v15 ==  v_axis15 
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
        axis *= cond_swap_BC(axis['B', 0, 8] == axis['C', 0, 8])
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
        l = [(1, i, (a[1,i] - 2) & 1) for i in (2, 3)]
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
        l = [(2, 3, (a[2, 3] - 1) & 1)]
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
    zi = axis.central_involution(guide) * axis.g_central
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

def beautify_axis(class_, g, verbose = 0, rand = 0, check = False):
    if REF_CLASS and class_ != REF_CLASS:
        # for debugging only
        return None

    axis = Axis(g)
    if verbose:
        print("\nOrbit" + class_)
        if not rand and mmv_scalprod(axis.v15, v_axis15) == 0:
            print("Axis may also be used in Baby monster")
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
            axis.display_sym_integral(e)
    return axis


#################################################################
# Main program
#################################################################


if __name__ == "__main__":
    from mmgroup.tests.axes.sample_axes import g_strings, g_classes
    for class_, g in zip(g_classes, g_strings):
        beautify_axis(class_, g, verbose = 1) 




