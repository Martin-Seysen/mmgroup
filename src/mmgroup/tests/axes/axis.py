r"""Investigate certain orbits of 2A axes of the monster group

The functions in this module are used to find the orbits of 2A axes
of the subgroup :math:`G_{x0}` (of structure 
:math:`2^{1+24}.\mbox{Co}_1`) of the monster group. According to
:cite:`Nor98` there are 12 such orbits.

"""

import sys
import os
import warnings
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

import mmgroup
from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2
from mmgroup.mm_crt_space import MMVectorCRT
from mmgroup.generators import mm_group_invert_word
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
from mmgroup.mm_op import mm_aux_mmv_size, mm_aux_index_mmv
from mmgroup.mm_op import mm_op_scalprod, mm_op_scalprod_ind




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

# Group element mapping v^+ to v^- and vice versa
G_MAP_STD_OPP = "x_200h"




G = MM0
G_SHORTEN = False
MM = None




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
    ERR =  "Cannot display a part of %s of a vector"
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
    r"""Return certain array of short vectors in Leech lattice mod 2

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
# Rebase an axis in the Monster
#################################################################




def rebase_axis(v15):
    try:
        from mmgroup.mm_reduce import mm_reduce_vector_vp
    except:
        W = "Cannot not rebase a monster axis"
        warnings.warn(W, UserWarning)
        return None
    v0 = np.zeros(1, dtype = np.uint32)
    v = v15.copy()
    w = MMV(15)(0)
    g = np.zeros(256, dtype = np.uint32)
    l_g = mm_reduce_vector_vp(v0, v.data, 1, g, w.data)
    assert 0 <= l_g < 128, hex(l_g)
    g0 = G('a', g[:l_g]) ** -1
    v_axis15 = Axis().v15
    if v_axis15 * g0 != v15:
        W = "Rebasing of a monster axis has failed"
        warnings.warn(W, UserWarning)
        return None
    return g0



#################################################################
# Sizes of G_x0 orbits on the axes in the Monster
#################################################################


AXIS_ORBIT_SIZES = {
'2A' : 196560,
'2B' : 11935123200,
'4A' : 1630347264000,
'4B' : 1466587938816000,
'4C' : 6599645724672000,
'6A' : 1896194506752000,
'6C' : 438020931059712000,
'8B' : 8601138282627072000,
'6F' : 1501786049347584000,
'10A' : 786389785840189440,
'10B' : 37845008443559116800,
'12C' : 48057153579122688000,
}



#################################################################
# Map scalar product of axes (mod 15) to product class
#################################################################

try:
   INDEX_BUFSIZE = mm_aux_mmv_size(15) // 2
except:
   INDEX_BUFSIZE = None

PROD_CLASSES = {
    8 : "1A",  1 : "2A",  14 : "3A",  4: "4A",  3 : "5A",
    10 : "6A",  0 : "2B",  2 : None
}

PROD_CLASSES_SQ = {1 : "4B", 2 : "3C"}



#################################################################
# Support for random axes
#################################################################

RandomAxis = None
RandomBabyAxis = None

def import_random_axis():
    global RandomAxis, RandomBabyAxis
    try:
        from mmgroup.tests.axes.random_axis import RandomAxis, RandomBabyAxis
    except:
        RandomAxis = RandomBabyAxis = None

#################################################################
#################################################################
# Class Axis for modelling a 2A axis in the Monster
#################################################################
#################################################################


class Axis:
    r"""Model a 2A axis in the representation of the Monster

    An instance of this class is an arbitrary 2A axis ``a`` of the
    Monster group with entries taken modulo 15. Attribute ``a.v15``
    of type |MMVector| is the vector in the representation of the
    Monster corrsponding to the axis. Apart from that vector we also
    maintain a property ``a.g`` of class |MM| such that the axis
    is ``a.v15`` equal to the product ``a.v15_start * a.g``. Here
    ``a.v15_start`` is the standard 2A axis :math:`v^+` defined in
    :cite:`Seysen22`. One may multiply an axis with an element of
    the Monster group. One may abbreviate ``a.v15[x]`` to ``a[x]``
    for any value ``x``.


    The constructor may of this class may be:

    * An element :math:`g` of the Monster. Then we constuct the
      axis :math:`v^+ \cdot g`.

    * The string ``'i'`` followed by a 2A involution :math:`g` in
      the Monster. Then we construct the axis corresponding
      to :math:`g`.

    * An instance of this class |MMVector| which is a 2A axis in
      representation of the Monster (mod 15).

    * An instance of class Axis. Then we construct a copy of that
      instance.
    """
    v15_start_name = V_AXIS 
    g_map_vstart_name = 1
    constant = False
    group = G
    ERR_CONST = "This axis is constant and may not be changed"
    ERR_MAP = "Could not map type %s object to a monster axis"
    ERR_INVOL = "Type %s object is not a 2A involution in the Monster"
    axis_class = "standard"
    default_profile_hashmode = 0
    def __init__(self, g = 1, invol = None):
        self.g0 = G()
        self.index_buf = None
        t = g
        if isinstance(g, Axis):
            self.g1 = G(g.g_map_vstart) * G(g.g)
        elif isinstance(g, MMVector):
            self.g1 = rebase_axis(g % 15)
        elif g == 'i':
            t = invol
            from mmgroup import MM
            g1 = MM(invol)
            i, h = g1.conjugate_involution()
            if i != 1:
                raise ValueError(ERR_INVOL % type(t))
            self.g1 = G(h**-1)
        elif g == 'r' and invol is None:
            # Not necessary, but much faster than the standard way
            if RandomAxis is None:
                import_random_axis()
            if RandomAxis is not None:
                ax = RandomAxis()
                self.g1 = G('a', ax.g('a'))
                self.v15 = ax.v15()
                return
            self.g1 = G('r')   # the standard way should still work
        elif g == 'a':
            # Undocumented function:  Axis('a', data) is equivalent
            # to  Axis(self.group('a', data)), but faster.
            self.g1 = G('a', invol)
        else:
            self.g1 = G(g)
        if self.g1 is None:
            raise ValueError(self.ERR_MAP % type(t))
        self.v15 = self.v_axis15 * self.g1
    @property
    def v_axis15(cls):
        r"""Return the fixed vector :math:`v^+` (mod 15)"""
        return MMV(15)(V_AXIS)
    @property
    def v15_start(cls):
        r"""Return the start axis vector of this class (mod 15)

        For an instance ``a`` of this class we have
        ``a.v15_start`` ==  ``a.v_axis15`` ==  :math:`v^+`. 
        """
        return MMV(15)(cls.v15_start_name)
    @property
    def g_central(cls):
        r"""Return central involution :math:`x_{-1}` of the group :math:`G_{x0}`"""
        return G(G_CENTRAL)
    @property
    def g_axis_start(cls):
        r"""Return the fixed 2A involution :math:`\beta^+`."""
        return G(G_AXIS)
    @property
    def g_map_vstart(cls):
        r"""Group element mapping :math:`v^+` to start axis

        For this class the neutral element of the group is returned.
        """
        return G(g_map_vstart_name)
    @classmethod
    def representatives(cls):
        r"""Map names :math:`G_{x0}` orbits of axes to their representatives

        The method returns an (ordered) dictionary that maps the names
        of the :math:`G_{x0}` orbits of the axes to the representatives
        of the orbits.
        """
        from mmgroup.tests.axes.get_sample_axes import get_sample_axes
        return get_sample_axes()
    @classmethod
    def orbit_sizes(cls):
        r"""Return sizes of the :math:`G_{x0}` orbits of axes

        The method returns dictionary that maps the names
        of the :math:`G_{x0}` orbits of the axes to their sizes.

        These orbit sizes are taken from :cite:`Nor98`.
        """
        return  AXIS_ORBIT_SIZES
    @property
    def g(self):
        r"""Group element mapping the start axis to the current axis 

        The function returns ``g`` with ``a.v15 = a.v15_start * a.g``.
        """
        return self.g0 * self.g1
    @property
    def g_axis(self):
        r"""Return the 2A involution corresponding to the 2A axis"""
        return self.g_axis_start ** self.g
    @property
    def group(self):
        r"""Return class implementing the Monster group acting on axes

        This is usually class |MM|, and need not be changed. One may
        change this by calling function ``set_axis_group``.
        """
        return G
    @property
    def norm_A_mod15(self):
        r"""Return norm of part A of the axis (modulo 15)."""
        return mm_op_norm_A(15, self.v15.data)
    def copy(self):
        r"""Return a (deep) copy of the axis"""
        axis = deepcopy(self)
        if axis.constant:
            del axis.constant
        return axis
    def __imul__(self, g):
        if self.constant:
            raise TypeError(self.ERR_CONST)
        g = G(g)
        self.g1 *= g
        self.v15 *= g
        self.index_buf = None
        return self
    def __mul__(self, g):
        return self.copy().__imul__(g)
    def __eq__(self, other):
        if not isinstance(other, Axis):
            ERR = "Cannot compare %s object with %s object"
            raise TypeError(ERR % (self.__class__, other.__class__))
        return self.v15 == other.v15
    def __getitem__(self, index):
        return self.v15.__getitem__(index)
    def rebase(self):
        r"""Rebase group element in axis

        For an instance ``a`` of class ``Axis`` we try to shorten the
        group element ``a.g`` with the full *mmgroup* functionality.
        Here we may also change ``a.g`` to a shorter element mapping
        the standard axis to the current axis.

        Class ``Axis`` will be used for proofs; so we always check
        the correctness of such an operation.

        We store ``a.g`` as a product ``a.g0 * a.g1``. Construction of
        ``a`` and multiplication of ``a`` with a group element affects
        ``a.g1`` only. This method puts ``a.g0 = a.g`` and ``a.g1 = 1``
        and reduces ``a.g0`` using the features of class |MM|.
        """
        if self.constant:
            raise TypeError(self.ERR_CONST)
        self.g0 *= self.g1
        self.g1 = G()
        g0 = rebase_axis(self.v15)
        if g0 is not None:
            self.g0 = g0
        return self
    def central_involution(self):
        r"""Return central involution of dihedral group

        The method returns the central involution of the dihedral
        group generated by ``a.g_axis`` and ``a.g_central``. For an
        instance ``a`` of this class, the group element ``a.g_axis``
        is the 2A involution corresponding to the axis; and
        ``a.g_central`` is the central involution in :math:`G_{x0}`.
        """
        c = self.g_axis * self.g_central
        _, h = c.half_order()
        return h
    def axis_type_info(self):
        r"""Return information about the :math:`G_{x0}` orbit of axis

        Let ``c = a.g_axis * a.g_central``. The method returns a triple
        ``(o, chi, s)``, where ``o`` is the (even) order of ``c``; and
        ``chi`` is the character of ``c`` in the representation of the
        Monster of dimension 196883. We cannot compute character ``chi``
        if the involution ``i = c**(o/2)`` is a 2A involution in the
        Monster; so we put ``chi = None`` in this case. The string ``s``   
        of the triple describes the class of involution ``i`` in
        :math:`G_{x0}` as in method ``conjugate_involution_G_x0`` of
        class ``mmgroup.Xsp2_Co1``.

        This information is sufficient to match the :math:`G_{x0}`
        orbit of an axis with one of the orbits given in
        :cite:`Nor98`.
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
        r"""Return array of short vectors in Leech lattice mod 2

        Let ``v = a.v15`` be the axis given by instance ``a`` of this
        class.  The function returns the list of all vectors ``v2`` in
        the Leech lattice mod 2 such that the absolute value of the
        entry of ``v`` corresponding a unit vector labelled by ``v2``
        is equal to ``value``. That list is returned as a numpy array
        with entries given in *Leech lattice encoding*.

        If ``radical`` is 1 then the function returns a basis of the
        radical of the space generated by the vectors described above.
        If ``radical`` is 2 then the function returns the whole
        radical as an array of vectors.
        """
        return find_short(self.v15, value, radical, verbose)
    def leech3matrix(self):
        r"""Return the 'A' part of the axis.

        The function returns the part of with tag 'A' of the axis
        as a matrix. That matrix a will be given in matrix mod 3
        encoding as documented in the header of
        file ``leech3matrix.c``.
        """
        a = np.zeros(3 * 24, dtype = np.uint64)
        mm_op_load_leech3matrix(self.v15.p, self.v15.data, a)
        return a
    def axis_type(self, e = 0):
        r"""Return the :math:`G_{x0}` orbit of the axis as a string

        ``a.axis_type(e)`` is equivalent to ``a.v15.axis_type(e)``.

        The list of possible axis types is obtained by calling
        ``Axis.representatives().keys()``
        """
        e %= 3
        v = self.v15.data
        if e:
            v = np.zeros(24*4, dtype = np.uint64)
            mm_op_t_A(15, self.v15.data, e, v)
        return axis_type(v)
    def in_space(self, space, *args):
        r"""Return axis as an element of a vector space ``space``

        The function returns the axis as an element of the class given
        by ``space(*args)``.
        """
        all_args = args + (self.v15_start_name,)
        return space(*all_args) * self.g
    def display_sym(self, part = 'A', mod = 0,
            diff = 0, ind = None, text = "", end = ""):
        r"""Display part of axis as a symmetric 24 times 24 matrix

        if ``part`` is 'A', 'B', or 'C' then the matrix corresponding to
        that part (modulo 15) is displayed. If ``part`` is an integer
        then the 'A' part of the matrix (multiplied with group element
        ``MM('t', part)``) is displayed. Parameter ``mod`` is 0 (default)
        or an optional odd modulus. If ``mod`` is odd then entries are
        displayed modulo ``mod``. If ``mod`` is negative then entries
        displayed modulo ``abs(mod)`` with absolute value less than
        ``mod / 2``. At present ``abs(mod)`` must be 3, 5, or 15.

        If ``diff`` is an instance of class ``Axis`` then ``diff``
        is subtracted from the axis. Parameter ``ind`` is an optional
        subset of the indices of entries of part ``A`` of the axis.

        ``text`` is an optional text to be displayed; ``end`` is an
        optional end mark.
        """
        if isinstance(part, str) and part.isdigit():
            part = int(part)
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
            v = self.in_space(MMVectorCRT, 20)
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
        r"""Compute invariants of a matrix related to axis part 'A'

        Let ``A`` be the symmetric 24 times 24 matrix corresponding
        to part 'A' of the axis, with entries taken modulo 3.
 
        Put ``M = A - d * E``, where  ``d`` is a integer and ``E`` is
        the unit matrix. We return a tuple ``(ker, isect, M_img, M_ker)``,
        where ``ker`` is the dimension of the kernel of ``M``, ``isect``
        is the dimension of the intersection of the kernel and the image
        of ``M``; and ``M_img`` and ``M_ker`` are 24 times 24 matrices.
        Row ``i`` of matrix ``M_img`` is the image of row ``i`` of
        ``M_ker`` under the symmetric matrix ``M``. The first ``isect``
        rows of ``M_img`` and ``M_ker`` are equal; they are a basis of
        the intersection of the kernel and the image of ``M``. The first
        ``24 - ker`` rows of ``M_img`` are a basis of the image of ``M``.
        The last ``ker`` rows of ``M_ker`` are a basis of the
        kernel of ``M``.
        """
        a = self.leech3matrix()
        leech3matrix_sub_diag(a, d, 0)
        leech3matrix_sub_diag(a, 2, 24)
        x = leech3matrix_kernel_image(a)
        leech3matrix_compress(a, a)
        img, ker, isect = x & 255, (x >> 8) & 255, x >> 16
        assert ker + img == 24
        return ker, isect, a[:24], a[24:48]
    def reduce_G_x0(self):
        r"""Compute a group element reducing the axis

        The function returns an element of the subgroup
        :math:`G_{x0}` of the Monster that maps the axis to the
        representative of its :math:`G_{x0}` orbit.
        """
        from mmgroup.tests.axes.reduce_axis import reduce_axis_G_x0
        return reduce_axis_G_x0(self)
    def profile_Nxyz(self, t = None, mode = None):
        r"""Return invariant of the axis under the group :math:`N_{xyz}`

        Let :math:`N_{xyz}` (of structure :math:`2^{2+11+22}.M_{24}`)
        be the unique subgroup of index 2 of :math:`N_{x0}`. The
        function computes a certain :math:`24 \times 24` matrix ``M``
        from the parts 'A', 'B', and 'C' of the axis modulo 3. Matrix
        ``M`` is invariant under the action of the normal subgroup
        :math:`2^{2+11+22}` of :math:`N_{xyz}`.

        The function also computes a matrix ``H`` corrsponding to the
        matrix ``M`` sorted in a  suitable way such that the sorted
        matrix ``H`` is invariant under :math:`N_{xyz}`. The function
        returns the triple ``(M, h, H)`` where the 63-bit integer ``h``
        is a hash value over ``H``.

        If parameter ``t`` is a tuple ``(e, f)`` then the same
        calculation is done with axis
        :math:`a \cdot \tau^e x_\delta^f` instead of the axis
        :math:`a` in this instance. Here :math:`\delta` may by any
        odd Golay cocode word. Note that
        :math:`\langle N_{xyz}, \tau^e, x_\delta^f
        \mid e=0,1,2; f=0,1 \rangle`
        is the maximal subgroup :math:`N_{0}` of the Monster of
        structure :math:`2^{2+11+22}.(S_3 \times M_{24})`.

        Parameter ``mode``  allows to specify a subgroup :math:`N'` of
        :math:`N_{xyz}`, such that the computed matrix ``M`` will be
        inviarant under  :math:`N'` only. Here :math:`N'` is the
        subgroup of :math:`N` centralizing a set
        :math:`\{x_\delta \mid \delta \in S\}` of 2A involutions
        pointwise, where :math:`S` is a set of Golay cocode words
        depending of parameter ``mode``. Legal values for ``mode`` are:

        .. table:: Subgroups :math:`N'` of :math:`N_{xyz}` supported
          :widths: 15 85

          ==== =============================================================
          mode :math:`N'` centralizes :math:`\{x_\delta \mid \delta \in S\}`
          ==== =============================================================
           0   :math:`S = \{ \}`, i.e. :math:`N' = N_{xyz}`
           1   :math:`S = \{ [2,3] \}`
           2   :math:`S = \{ [1,2], [2,3] \}`
          ==== =============================================================
        """
        from mmgroup.mm_reduce import mm_profile_mod3_load
        from mmgroup.mm_reduce import mm_profile_mod3_hash
        if isinstance(t, tuple):
            t = (t[0] % 3) * 2 + (t[1] & 1)
        else:
            t = 0
        if mode is None:
            mode = self.default_profile_hashmode
        a = np.zeros(72, dtype = np.uint64)
        b = np.zeros(2*576, dtype = np.uint16)
        assert  mm_profile_mod3_load(15, self.v15.data, a, t) == 0
        h = int(mm_profile_mod3_hash(a, b, mode))
        assert h >= 0
        b = b.reshape((2,24,24))
        #print(b)
        return b[0], h, b[1]
    def scalprod15(self, axis, sparse = False):
        r"""Scalar product of two axes (modulo 15)

        This is an auxiliary method for method ``product_class``. It
        returns the scalar product of this axis and another axis given
        by parameter ``axis`` (modulo 15). Parameter ``sparse`` is as
        in method ``product_class``.
        """
        v1 = self.v15.data
        v2 = axis.v15.data
        if sparse and self.index_buf is None:
            buf = np.zeros(INDEX_BUFSIZE, dtype = np.uint16)
            n1 = mm_aux_index_mmv(15, v1, buf, INDEX_BUFSIZE)
            self.index_buf = buf[:n1].copy() if n1 > 0 else False
        if isinstance(self.index_buf, np.ndarray):
            res = mm_op_scalprod_ind(15, v1, v2, self.index_buf)
        else:
            res = mm_op_scalprod(15, v1, v2)
        assert res >= 0
        return res
    def product_class(self, axis, sparse = False):
        r"""Conjugation class of product of this axis with another axis

        The method returns the conjugation class of the product of
        the involution corresponding to this axis with the involution
        corresponding to another axis given by parameter ``axis``.
        This conjugation class (in the Monster) is one of the strings
        '1A', '2A', '2B', '3A', '3C', '4A', '4B', '5A', or '6A'. It
        is computed using scalar products of axes, see e.g.
        :cite:`Con85`, Fig. 3.

        If parameter ``sparse`` is ``True`` then some optimizations
        are made for the case that this axis is a sparse vector.
        """
        cl = PROD_CLASSES[self.scalprod15(axis, sparse)]
        if cl is None:
            v = self * axis.g_axis
            cl = PROD_CLASSES_SQ[self.scalprod15(v)]
        return cl
    def _fast_g_mmdata(self, e = 1):
        r"""Equivalent to (self.g**e).mmdata for e = +-1

        This is faster than the standard procedure and makes
        no attempt to reduce the returnd data.
        """
        assert abs(e) == 1
        data = np.concatenate((self.g0.mmdata, self.g1.mmdata))
        if e == -1:
            mm_group_invert_word(data, len(data))
        return data


#################################################################
#################################################################
# Auxiliary function for class BabyAxis
#################################################################
#################################################################


#################################################################
# Constructive membership test of Monster element in Baby Monster
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

def in_H(g):
    try:
        mm, lmm = g.mmdata, len(g.mmdata)
        for x in [V_PLUS, 0x1000000]:
            if gen_leech2_op_word(x, mm, lmm) != x:
                return None
        g1 = G(Xsp2_Co1(g))
        return in_Baby_direct(g1)
    except:
        return None
 
def in_Baby(g):
    """Constructive membership test in Baby Monster ``2.B``

    Given an element ``g`` of the Monster we try to find a word ``g1``
    of generators of the Baby Monster with ``g == g1``. If  we can
    find such a ``g1`` then we return ``g1`` as an instance of
    class ``Axis.group``. Otherwise we return ``None``.
    """
    if in_Baby_direct(g) is not None:
        return g
    g1 = in_H(g)
    if g1 is not None:
        return g1
    try:
        from mmgroup import MM
        v_axis15 = Axis().v15
        v_axis_opp15 = BabyAxis().v15
        g1 = G(MM(g).reduce())
        if v_axis15 * g1 == v_axis15:
            if v_axis_opp15 * g == v_axis_opp15 * g1:
                return in_Baby_direct(g)
        return None
    except:
        return None



#################################################################
# Rebase an axis in the Baby Monster
#################################################################




def rebase_baby_axis(v15):
    try:
        from mmgroup.mm_reduce import mm_reduce_vector_shortcut
        from mmgroup.mm_reduce import mm_reduce_vector_vm
        from mmgroup import mmv_scalprod
    except:
        W = "Cannot not rebase a baby monster axis"
        warnings.warn(W, UserWarning)
        return None
    v0 = np.zeros(1, dtype = np.uint32)
    v = v15.copy()
    w = MMV(15)(0)
    g = np.zeros(256, dtype = np.uint32)
    l_g = mm_reduce_vector_shortcut(1, 1, V_PLUS, g)
    assert 0 <= l_g < 128
    l_g = mm_reduce_vector_vm(v0, v.data, g, w.data)
    assert 0 <= l_g < 128, hex(l_g)
    g0 = G('a', g[:l_g]) ** -1
    v_axis15 = Axis().v15
    v_axis_opp15 = BabyAxis().v15
    if (v_axis15 * g0 != v_axis15 or v_axis_opp15 * g0 != v15 or
        mmv_scalprod(v15, v_axis15) != 0):
        W = "Rebasing of a baby monster axis has failed"
        warnings.warn(W, UserWarning)
        return None
    return g0

#################################################################
# Sizes of G_x0 \cap 2.B orbits on the axes in the Baby Monster
#################################################################

BABY_AXIS_ORBIT_SIZES = {
'2A1' : 1,
'2A0' : 93150,
'2B1' : 7286400,
'2B0' : 262310400,
'4A1' : 4196966400,
'4B1' : 470060236800,
'4C1' : 537211699200,
'6A1' : 9646899200,
'6C1' : 6685301145600,
'10A1' : 4000762036224,
}

#################################################################
#################################################################
# Class BabyAxis for modelling a 2A axis in the Baby monster
#################################################################
#################################################################




class BabyAxis(Axis):
    r"""Model a 2A axis in the representation of the Baby Monster

    The constructor and the methods of this class are as in the base
    class ``Axis`` of this class. But we accept axes that are
    orthogonal to axis :math:`v^+` only. The centralizer of that
    axis is a subgroup of structure :math:`2.B` of the Monster.

    For an instance ``a`` of this class we have
    ``a.v15 = a.v15_start * a.g``, where ``a.v15_start`` is the
    (fixed) axis :math:`v^-` opposite to the start axis :math:`v^+`
    of class ``Axis``.
    """
    v15_start_name = V_AXIS_OPP
    g_map_vstart_name = 'x_0x200'
    constant = False
    ERR_BABY_G = "Type %s object is not in the baby monster group"
    ERR_BABY = "Could not map type %s object to a baby monster axis"
    axis_class = "baby"
    default_profile_hashmode = 1
    def __init__(self, g = 1, invol = None):
        self.g0 = G()
        self.index_buf = None
        t = g
        if isinstance(g, BabyAxis):
            self.g1 = g.g
        elif isinstance(g, Axis):
            self.g1 = rebase_baby_axis(g.v15)
        elif g == 'i':
            v15 = Axis('i', invol).v15
            self.g1 = rebase_baby_axis(v15)
            t = invol
        elif isinstance(g, MMVector):
            self.g1 = rebase_baby_axis(g % 15)
        elif g == 'r' and invol is None:
            # Not necessary, but much faster than the standard way
            if RandomBabyAxis is None:
                import_random_axis()
            if RandomBabyAxis is not None:
                ax = RandomBabyAxis()
                self.g1 = G('a', ax.g('a'))
                self.v15 = ax.v15()
                return
            self.g1 = G('r', 'B') # the standard way should still work
        else:
            self.g1 = in_Baby(G(g))
        if self.g1 is None:
            raise ValueError(self.ERR_BABY % type(t))
        self.v15 = self.v15_start * self.g
    @property
    def g_map_vstart(cls):
        r"""Group element mapping :math:`v^+` to the start axis

        For this class the start axis is :math:`v^-`.
        """
        return G(cls.g_map_vstart_name)
    @classmethod
    def representatives(cls):
        r"""Map names :math:`H` orbits of axes to their representatives

        The method returns an (ordered) dictionary that maps the names
        of the :math:`H` orbits of the axes to the representatives
        of the orbits. Here :math:`H = G_{x0} \cap H^+`, where
        :math:`H^+` (of structure :math:`2.B`) is the centralizer
        of the standard axis :math:`v^+`.
        """
        from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
        return get_baby_sample_axes()
    @classmethod
    def orbit_sizes(cls):
        r"""Return sizes of the :math:`H` orbits of axes

        The method returns dictionary that maps the names
        of the :math:`H` orbits of the axes to their sizes.
        Here :math:`H = G_{x0} \cap H^+`, where :math:`H^+` 
        (of structure :math:`2.B`) is the centralizer
        of the standard axis :math:`v^+`.

        These orbit sizes are taken from :cite:`Mue08`.
        """
        return BABY_AXIS_ORBIT_SIZES
    @property
    def g_axis_start(self):
        r"""Return the fixed 2A involution :math:`\beta^-` """
        return G(G_AXIS_OPP)
    def __imul__(self, g):
        if self.constant:
            raise TypeError(self.ERR_CONST)
        g = in_Baby(G(g))
        if g is None:
            raise ValueError(self.ERR_BABY_G % type(g))
        self.g1 *= g
        self.v15 *= g
        return self
    def __mul__(self, g):
        return self.copy().__imul__(g)
    def rebase(self):
        if self.constant:
            raise TypeError(self.ERR_CONST)
        self.g0 *= self.g1
        self.g1 = G()
        g0 = rebase_baby_axis(self.v15)
        if g0 is not None:
            self.g0 = g0
        return self
    def fixed_value(self, part = 'A'):
        if part == 'A':
            a = self['A']
            return (int(a[2,2]) + int(a[3,3]) - 2 * int(a[2,3])) % 15
        elif part == 'B':
            return int(self['B', 2, 3])
        else:
            E = "Illegal argument for method fix_value of class BabyAxis"
            raise ValueError(E)
    def axis_type(self, e = 0):
        r"""Return the :math:`H` orbit of the axis as a string


        Here :math:`H = G_{x0} \cap H^+`, where  :math:`H^+` (of
        structure :math:`2.B`) is the centralizer of the standard
        axis :math:`v^+`.

        ``a.axis_type(e)`` is equivalent to ``a.v15.axis_type(e)``.

        The list of possible axis types is obtained by calling
        ``BabyAxis.representatives().keys()``
        """
        at =  super(BabyAxis, self).axis_type(e)
        if e:
            vt = self.v15.copy()
            mm_op_t_A(15, self.v15.data, e % 3, vt.data)
        else:
            vt = self.v15
        a = vt['A']
        asub = (int(a[2,2]) + int(a[2,3]) - 2 * int(a[2,3])) % 15
        return at + str(int(asub != 0)) 
    def axis_type_info_beta(self):
        r"""Return information about a certain orbit of the axis

        Let ``c = a.g_axis * a.g_central``; and let ``o`` be the (even)
        order of ``c``. Let ``i = c**(o/2)``. Then the involution 
        ``i1 = i * Axis.g_axis_start()`` is in  :math:`G_{x0}`.
        We return the string ``s``  describing the class of
        involution  ``i`` in :math:`G_{x0}`, as in method
        ``conjugate_involution_G_x0`` of class ``mmgroup.Xsp2_Co1``.
        """
        c = self.g_axis * self.g_central
        o, i = c.half_order()
        assert o & 1 == 0
        i *= G(G_AXIS)
        s, _ = i.conjugate_involution_G_x0()
        return s
    def reduce_G_x0(self):
        r"""Compute a group element reducing the axis

        Let :math:`H` be the the intersection of the centralizer of
        the standard 2A axis :math:`v^+` (of structure :math:`2.B`)
        with the subgroup :math:`G_{x0}` of the Monster. 
        The function returns an element of :math:`H` that maps the
        axis to the representative of its :math:`H` orbit.
        """
        from mmgroup.tests.axes.reduce_baby_axis import reduce_baby_axis_G_x0
        return reduce_baby_axis_G_x0(self)





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

    Classes ``Axis`` and ``BabyAxis`` are used for proofs. Therefore
    it may be useful to compute in the Monster using instances of class
    ``MM0`` instead of class |MM|. Class ``MM0`` implements operations
    of the Monster defined in :cite:`Seysen20` only. Note that in class
    ``MM0`` we cannot shorten the word representing an element of the
    Monster. One may use method ``rebase`` of class ``Axis`` for
    shortening the word representing a group element in the axis.

    If users need this function they should call it once and
    for all before dealing with any 2A axes.
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


if not MM_INITIALIZED:
    set_axis_group()



