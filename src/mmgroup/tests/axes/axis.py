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

from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup.generators import gen_leech3to2
from mmgroup.generators import gen_leech3to2_short
from mmgroup.generators import gen_leech2_op_atom, gen_leech2_op_word
from mmgroup.clifford12 import leech3matrix_kernel_image
from mmgroup.clifford12 import leech3matrix_compress
from mmgroup.clifford12 import leech3matrix_sub_diag
from mmgroup.clifford12 import leech2_matrix_basis
from mmgroup.clifford12 import leech2_matrix_radical
from mmgroup.clifford12 import leech2_matrix_expand
from mmgroup.mm_op import mm_op_load_leech3matrix
from mmgroup.mm_op import mm_op_eval_X_find_abs
from mmgroup.mm_op import mm_op_norm_A
from mmgroup.mm_op import mm_op_eval_A
from mmgroup.mm_op import mm_op_t_A
from mmgroup.mm_op import mm_op_eval_A_rank_mod3




########################################################################
########################################################################
# Start axis vectors and involutions in the monster group
########################################################################
########################################################################

# We give the standard axes and involutions as strings. So we can use
# it in any suitable constructor of the for the monster or its rep.

# Standard axes v^+ of 2A involution x_\beta, \beta = Cocode([2,3])
V_AXIS = "A_2_2 - A_3_2 + A_3_3 - 2*B_3_2" 
# Opposite axis v^- of of 2A involution x_{-1} * x_\beta  
V_AXIS_OPP = "A_2_2 - A_3_2 + A_3_3 + 2*B_3_2" 

# 2A involution x_\beta corresponding to standard axis v^+
G_AXIS = "d_200h"
# 2A involution x_{-1} x_\beta corresponding to opposite axis v^-
G_AXIS_OPP = "x_1000h * d_200h"

# Central involution in the subgroup G_x0 of the monster
G_CENTRAL = "x_1000h"

# Group element mapping v^+ to v^-
G_MAP_STD_OPP = "x_200h"



V15 = MMV(15)
v_axis = MMVectorCRT(16, V_AXIS)
v_axis_opp = MMVectorCRT(16, V_AXIS_OPP)
v_axis15 = V15(V_AXIS)
v_axis_opp15 = V15(V_AXIS_OPP)


G = MM0
G_SHORTEN = False
MM = None

# The following constants will be set by function set_axis_group

# The central involution in the subgroup ``G_x0``-
# g_central = G(G_CENTRAL)  

# The standard 2A element in the monster group
# g_axis = G(G_AXIS)

# The opposite standard 2A element in the monster group
# g_axis_opp = G(G_AXIS_OPP)

# Element mapping the standard axis to its opposite axis
# g_map_std_opp = G(G_MAP_STD_OPP)




#################################################################
#################################################################
# Auxiliary function for class Axis
#################################################################
#################################################################



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


def common_factor_v2(A):
    z = reduce(__or__, np.array(A.ravel(), dtype = np.uint32), 0)
    if z & 255 == 0:
        return 256
    if z & 15 == 0:
        return 16
    return 1

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
# Determine the  G_x0 orbit of an axis
#################################################################


_NORM_TYPES = {2:"8B", 3:"4C", 5:"6C", 10:"12C", 13:"4B"}
_RANK_TYPES = {(8, 8): "2B", (8, 24): "10B",
               (14, 8): "6F", (14, 23): "4A"}
def axis_type(v):
    norm = mm_op_norm_A(15, v)
    if norm in _NORM_TYPES:
        return _NORM_TYPES[norm]
    if norm == 4:
        r =  mm_op_eval_A_rank_mod3(15, v, 2);
        rank = r >> 48;
        if rank == 2:
            return "10A"
        if rank == 23:
            v2 = gen_leech3to2_short(r & 0xffffffffffff) & 0xffffff
            val_A = mm_op_eval_A(15, v, v2);
            if val_A == 4:
                return "2A"
            if val_A == 7:
                return "6A"
    rank = mm_op_eval_A_rank_mod3(15, v, 0) >> 48
    try:
        return _RANK_TYPES[(norm, rank)]
    except:
        ERR = "Vector is not a 2A axis"
        raise valueError(ERR)



#################################################################
# Analyse the part of an axis that encodes short vectors 
#################################################################


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
#################################################################
# Class Axis for modelling a 2A axis in the Monster
#################################################################
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
    def axis_type_info(self):
        """Return information about the math:`G_{x0}` orbit of axis

        Let ``c = a.g_axis * a.g_central``. The method returns a triple
        ``(o, chi, s)``, where ``o`` is the (even) order of ``c``; and
        ``chi`` is the character of ``c`` in the rep of the Monster of
        dimension 196883. We cannot compute character ``chi`` if the
        involution ``i = c**(o/2)`` is a 2A involution in the Monster;
        so we put ``chi = None`` in this case. The string ``s`` of the
        triple describes the class of involution ``i`` in :math:`G_{x0}`
        as in method ``conjugate_involution_G_x0`` of class
        ``mmgroup.Xsp2_Co1``.

        This information is sufficient to match the :math:`G_{x0}`
        orbit of an axis with one of the orbits given in |Nor98|.
        """
        c = self.g_axis * self.g_central
        o, i = c.half_order()
        assert o & 1 == 0
        s, _ = i.conjugate_involution_G_x0()
        if i.chi_G_x0()[0] != 275:
            return o, None, s # Then  i = c**(o/2)  is a 2A involution
        # Here ``i`` is a 2B involution
        itype, m = i.conjugate_involution()
        assert itype == 2
        chi = (c**m).chi_G_x0()[0]
        return o, chi, s
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
        e %= 3
        v = self.v15.data
        if e:
            v = np.zeros(24*4, dtype = np.uint64)
            mm_op_t_A(15, self.v15.data, e, v)
        return axis_type(v)
    def axis_in_space(self, space, *args):
        """Return axis as an element of a ``space``

        The function returns the axis as an element of the space given
        by ``space(*args)``.
        """
        all_args = args + (self.v15_start_name,)
        return space(*all_args) * self.g
    def display_sym(self, part = 'A', mod = 0,
            diff = 0, ind = None, text = "", end = ""):
        """Display part of axis as a symmetric 24 times 24 matrix

        if ``part`` is 'A', 'B', or 'C' then the matrix corresponding to
        that part (modulo 15) is displayed. If ``part`` is an integer
        then the 'A' part of the matrix (multiplied with group element
        ``MM('t', part)``) is displayed. Parmeter ``mod`` is 0 (default)
        or an optional odd modulus. If ``mod`` is odd then entries are
        displayed modulo ``mod``. If ``mod`` is negative then entries
        displayed modulo ``abs(mod)`` with absolute value less than
        ``mod / 2``. At present ``abs(mod)`` must be 3, 5, or 15.

        If ``diff`` is is an instance of class ``Axis`` then ``diff``
        is subtracted from the axis. Parameter ``ind`` is an optional
        subset of the indices of entries of part ``A`` of the axis.

        ``text`` is an optional text to be displayed; ``end`` is an
        optional end mark.
        """
        if isinstance(part, Integral):
            part %= 3
        if not text:
            if not part:
                part = 'A'
            if part in [1, 2]:
                text = "Axis part 'A' * t**%d" % part
            elif part in "ABC":
                text = "Axis part '%s'" % part
        signed, mod = mod < 0, abs(mod)
        if mod:
            if mod not in [3, 5, 15]:
                ERR = "Display modulus must be 3, 5, or 15"
                raise ValueError(ERR)
            a =  np.array(sym_part(self.v15, part), dtype = np.int32)
            if diff:
                a = (a - sym_part(diff.v15, part))
            a = np.array(a % mod)
            assert a.dtype == np.int32
            if signed:
                a = np.choose(a <= mod//2, [a - mod, a])
            text += " (mod %d)" % mod
        else:
            if diff:
                ERR = "Parameter 'diff' is illegal if mod is 0"
                raise ValueError(ERR)
            v = self.axis_in_space(MMVectorCRT, 20)
            a = 256 * sym_part(v, part)
            f = common_factor_v2(a)
            a = a // f
            text += " (multiplied by %d)" % (256 // f)
        if ind is not None:
            a = np.array([[a[i, j] for j in ind] for i in ind])
        print(text)
        display_A(a)
        print(end)
    def kernel_A(self, d):
        """Compute invariants of a matrix related to axis part 'A'

        Let 'A' be the symmetric 24 times 24 matrix corresponding to
        part 'A' of the axis, with entries taken modulo 3.
 
        Put 'M' = 'A' - 'd' * 'E', where  'd' is a integer and 'E' is
        the unit matrix. We return a tuple '(ker, isect, M_img, M_ker)',
        where 'ker' is the dimension of the kernel of 'M', 'isect' is
        the dimension of the intersection of the kernel and the image
        of 'M'; and 'M_img' and 'M_ker' are 24 times 24 matrices. Row
        'i' of matrix 'M_img' is the image of row 'i' of 'M_ker' under
        the symmetric matrix 'M'. The first 'isect' rows of 'M_img'
        and 'M_ker' are equal; they are a basis of the intersection of
        the kernel and the image of 'M'. The first '24 - ker' rows of
        'M_img' are a basis of the image of 'M' The last 'ker' rows
        of 'M_ker' are a basis of the kernel of 'M'.
        """
        a = self.leech3matrix()
        leech3matrix_sub_diag(a, d, 0)
        leech3matrix_sub_diag(a, 2, 24)
        x = leech3matrix_kernel_image(a)
        leech3matrix_compress(a, a)
        img, ker, isect = x & 255, (x >> 8) & 255, x >> 16
        assert ker + img == 24
        return ker, isect, a[:24], a[24:48]




#################################################################
#################################################################
# Auxiliary function for class BabyAxis
#################################################################
#################################################################





V_PLUS = 0x200

def in_Baby_direct(g):
    for x in g.mmdata:
        mask = (x >> 28) & 7
        if mask == 2:  # case x_pi
            if gen_leech2_op_atom(V_PLUS, x) != V_PLUS:
                return None
        if mask in [3,4] and x & 0x200:
            return None
        if mask == 7:
            return None
    return g

H_AXES = [0x200, 0x1000000]
def in_H(g):
    mm, lmm = g.mmdata, len(g.mmdata)
    if any([gen_leech2_op_word(x, mm, lmm) != x for x in H_AXES]):
         return None
    g1 = G(Xsp2_Co1(g))
    assert in_Baby_direct(g1)
    return g1
 
def in_Baby(g):
    if in_Baby_direct(g) is not None:
        return g
    return in_H(g)




#################################################################
#################################################################
# Class BabyAxis for modelling a 2A axis in the Baby monster
#################################################################
#################################################################




class BabyAxis(Axis):
    v15_start = v_axis_opp15
    constant = False
    ERR_BABY = "Cannot map element %s to subgroup 2.B of Monster"
    """Models an axis modulo 15"""
    def __init__(self, g = 1):
        if isinstance(g, Axis):
            self.g1 = in_Baby(Axis.g)
        else:
            self.g1 = in_Baby(G(g))
        if self.g1 is None:
            raise ValueError(self.ERR_BABY % self.g1)
        self.g0 = G()
        self.v15 = self.v15_start * self.g
    @property
    def g_axis_start(self):
        """Return the fixed 2A axis ``v^-``"""
        return G(G_AXIS_OPP)
    def __imul__(self, g):
        if self.constant:
            raise TypeError(self.ERR1)
        g = in_Baby(G(g))
        if g is None:
            raise ValueError(self.ERR_BABY % g)
        self.g1 *= g
        self.v15 *= g
        return self
    def __mul__(self, g):
        return self.copy().__imul__(g)
    #def __getitem__(self, index):
    #    return self.v15.__getitem__(index)
    def rebase(self):
        if self.constant:
            raise TypeError(self.ERR1)
        self.g0 *= self.g1
        self.g1 = G()
        if v_axis * self.g0 != v_axis:
             ERR = "Involution of 2A axis is not in Baby Monster"
             raise ValueError(ERR)
        try:
            #print("Rebasing axis in Baby Monster")
            from mmgroup.mm_reduce import mm_reduce_vector_shortcut
            from mmgroup.mm_reduce import mm_reduce_vector_vm
            v0 = np.zeros(1, dtype = np.uint32)
            v = self.v15.copy()
            w = V15(0)
            g = np.zeros(128, dtype = np.uint32)
            l_g = mm_reduce_vector_shortcut(1, 1, V_PLUS, g)
            assert 0 <= l_g < 128
            l_g = mm_reduce_vector_vm(v0, v.data, g, w.data)
            assert 0 <= l_g < 128, hex(l_g)
            g0 = G('a', g[:l_g]) ** -1
            assert self.v15_start * g0 == self.v15
            self.g0 = g0
        except:
            import warnings
            W = "Reducing a baby monster axis with mmgroup has failed"
            warnings.warn(W, UserWarning)
            raise
        return self
    def fixed_value(self, part = 'A'):
        if part == 'A':
            a = self['A']
            return int(a[2,2] + a[2,3] - 2 * a[2,3]) % 15
        elif part == 'B':
            return int(self['B', 2, 3])
        else:
            E = "Illegal argument for method fix_value of class BabyAxis"
            raise ValueError(E)
    def axis_type(self, e = 0):
        at =  super(BabyAxis, self).axis_type(e)
        if e:
            vt = self.v15.copy()
            mm_op_t_A(15, self.v15.data, e % 3, vt.data)
        else:
            vt = self.v15
        a = vt['A']
        asub = int(a[2,2] + a[2,3] - 2 * a[2,3]) % 15
        return at + str(int(asub != 0)) 




#################################################################
#################################################################
# initialization of the module
#################################################################
#################################################################

MM_INITIALIZED = False


def set_axis_group(group = None, shorten = True):
    """Set group acting on instances of class ``Axis``.

    The ordinary user almost never has a resaon to call this
    function. Here parameter ``group`` should be ``None``
    (default), or  one of classes ``MM`` or ``MM0`` defined 
    im module ``mmgroup`. By default we use class ``MM`` if
    present and class ``MM0`` otherwise.

    More details about use cases for class ``MM0`` are yet
    to be documented!

    If users need this function they should call it once and
    for all before dealing with any 2A axes.
    """ 
    global G, G_SHORTEN, MM, MM_INITIALIZED, Axis
    global g_central, g_axis, g_axis_opp, g_map_std_opp
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

    # The central involution in the subgroup ``G_x0``-
    g_central = G(G_CENTRAL)  

    # The standard 2A element in the monster group
    g_axis = G(G_AXIS)

    # The opposite standard 2A element in the monster group
    g_axis_opp = G(G_AXIS_OPP)

    # Element mapping the standard axis to its opposite axis
    g_map_std_opp = G(G_MAP_STD_OPP)

    assert g_axis ** g_map_std_opp == g_axis_opp
    assert v_axis * g_map_std_opp == v_axis_opp
 
    MM_INITIALIZED = True


if not MM_INITIALIZED:
    set_axis_group()



