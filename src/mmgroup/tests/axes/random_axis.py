"""We implement Random 2A involutions in the Monster and their axes.

Class ``Axis`` in module ``mmgroup.axes`` models a 2A axis for the 
purpose of inspecting the axis and the 2A involutions corresponding to
them. There are applications of axes where one has to generate random
axes with utmost speed. E.g. for computing normalizers of small
subgroup of the Monster one has to generate a large number of axes and
do computionas with them, see :cite:`BPW16` for background.

The class ``RandomAxis`` implemented here allows very fast generation
of a random axis and of the corresponding 2A involution.

Another application of the functions in this module is fast
generation of truly uniform distributed elements of the Monster.

More details will be documented in a future version.
 
"""

# If import fails then we want a quick failure
import mmgroup
from mmgroup.tests.axes import sample_axes
from mmgroup.tests.axes import baby_sample_axes


from random import randint, uniform

import numpy as np

import mmgroup
from mmgroup import MM0, MM, MMV, MMVector, Cocode, XLeech2
from mmgroup import MAT24_ORDER, Xsp2_Co1, characteristics
try:
   from mmgroup import mat24
except (ImportError, ModuleNotFoundError):
   from mmgroup.dev.mat24.mat24_ref import Mat24
   mat24 = Mat24
from mmgroup.generators import gen_rng_uniform_to_dist
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import mm_group_invert_word
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2 
from mmgroup.mm_op import mm_vector, mm_op_word
from mmgroup.mm_reduce import gt_word_shorten 
from mmgroup.mm_reduce import mm_compress_expand_type4 
from mmgroup.tests.axes.axis import Axis, BabyAxis
from mmgroup.tests.axes.axis import G_AXIS, G_AXIS_OPP
from mmgroup.tests.axes.axis import V_AXIS, V_AXIS_OPP

characteristics()
VV = MMV(15)

#############################################################
# Random elements in G_x0 and in G_x0 \cap 2.B
#############################################################

MAX_LEN_G_x0 = 10
MAX_NUM_G_x0 = MAT24_ORDER << (25 + 11 + 23)

def rand_G_x0(num = None):
    r"""Generate Random element of G_x0 as array of 32-bit integers
    """
    a = np.zeros(MAX_LEN_G_x0, dtype = np.uint32)
    if num is None:
        y = randint(0, 0x7ff)
        x = randint(0, 0x1fff)
        d = randint(0, 0xfff)
        p = randint(0, MAT24_ORDER - 1) 
        c = -1
        while c < 0:
            c = mm_compress_expand_type4(randint(0, 0x7fffff))
        assert gen_leech2_type(c) == 4
        num = ((((c * MAT24_ORDER + p) * 0x1000) + d) * 0x2000) + x
        num = num * 0x800 + y
    else:
        num = int(num)
        num1, y  = divmod(num, 0x800)
        num1, x  = divmod(num1, 0x2000)
        num1, d  = divmod(num1, 0x1000)
        c, p = divmod(num1, MAT24_ORDER)
        assert 0 <= c < 0x1000000
    a[0] =  0x40000000 + y  # tag y
    a[1] =  0x30000000 + x  # tag x
    a[2] =  0x10000000 + d  # tag d
    a[3] =  0x20000000 + p  # tag p
    # Append a word in G_x0 that maps Omega to c
    len_a = gen_leech2_reduce_type4(c, a[4:])
    assert 0 <= len_a <= MAX_LEN_G_x0  - 4    
    mm_group_invert_word(a[4:], len_a)
    return a[:4 + len_a], num


MAT22_DBL_ORDER = MAT24_ORDER // (12 * 23)
MAX_NUM_H = MAT22_DBL_ORDER << 36

BETA = 0x200
MODE_MAT22 = 1  # mode for generating random element of Mat_22:2

def rand_G_x0_baby(num = None):
    r"""Random element of G_x0 \cap 2.B as array of 32-bit integers
    """
    a = np.zeros(MAX_LEN_G_x0, dtype = np.uint32)
    if num is None:
        y = randint(0, 0x7ff)
        x = randint(0, 0x1fff)
        d = randint(0, 0xfff)
        p = randint(0, MAT22_DBL_ORDER - 1) 
        c, v4 = 0, 0
        while gen_leech2_type(v4) != 4:
             # generate short vector v2 in Leech lattice mod 2
             c = randint(300, 98579)  
             vs = mm_aux_index_extern_to_sparse(c)
             v2 = mm_aux_index_sparse_to_leech2(vs) 
             # Check if v2 is orthogonal to BETA in the real Leech lattice
             v4 = v2 ^ BETA
        num = ((((c * MAT22_DBL_ORDER + p) * 0x1000) + d) * 0x2000) + x
        num = num * 0x800 + y
    else:
        num = int(num)
        num1, y  = divmod(num, 0x800)
        num1, x  = divmod(num1, 0x2000)
        num1, d  = divmod(num1, 0x1000)
        c, p = divmod(num1, MAT22_DBL_ORDER)
        assert 300 <= c <= 98579
        vs = mm_aux_index_extern_to_sparse(c)
        v2 = mm_aux_index_sparse_to_leech2(vs) 
        v4 = v2 ^ BETA
        assert gen_leech2_type(v4) == 4
    pi = mat24.m24num_rand_local(MODE_MAT22, p)
    a[0] =  0x40000000 + (y & ~BETA)  # tag y
    a[1] =  0x30000000 + (x & ~BETA)  # tag x
    a[2] =  0x10000000 + d  # tag d
    a[3] =  0x20000000 + pi  # tag p
    # Append a word in G_x0 that maps Omega to c
    len_a = gen_leech2_reduce_type4(v4, a[4:])
    assert 0 <= len_a <= MAX_LEN_G_x0  - 4, len_a   
    mm_group_invert_word(a[4:], len_a)
    return a[:4 + len_a], num






#############################################################
# Auxiliary functions
#############################################################



def chop_axis_representative(h):
    """Chop representative h of a G_x0 orbit of an axis

    Here h is usually an element of M such that v = v_0 * h
    represents the G_x0 orbit of an axis, where v_0 is the
    standard axis in the Monster. To generate a random axis,
    we will multiply v (or h) with a random element of G_x0.

    In mmgroup, a representative v (or h) of a G_x0 orbit of an
    axis as above has been selected because of its beauty, not
    because of its shortness. To obtain a random axis in its G_x0
    orbit, we will right multiply h with a random element of G_x0.
    Therefore we may delete all trailing generators of h that
    are in G_x0. This function deletes all these trailing
    generators from  h and returns the result.
    """
    cls = h.__class__
    tuples = h.as_tuples()
    while len(tuples) and tuples[-1][0] != 't':
        del tuples[-1]
    return cls(tuples)
    

def aggressive_reduce_word(g):
    """Shorten a word representing an element g of the Monster

    This function is usually applied to a word  g = h**-1 * i_0 * h,
    where h represents a G_x0 orbit of an axis as in function
    chop_axis_representative(); and i_0 is the 2A involution in
    the Monster corresponding to the standard axis v_0.

    We shorten the word representing  g  by applying the known
    relations in the subgroups G_x0 and N_0 of the Monster. In most
    cases the leads to a shorter word representing g than the
    standard reduction process in mmgroup.
    """
    cls = g.__class__
    gbuf = g.mmdata
    buf = np.zeros(256, dtype = np.uint32)
    l_buf = gt_word_shorten(gbuf, len(gbuf), buf, len(buf), mode=3)
    assert l_buf >= 0
    result = cls('a', buf[:l_buf])
    assert MM(result) == MM(g)
    return result

#############################################################
# Initializng the module
#############################################################



MAX_LEN_G = 40
MAX_LEN_INVOLUTION = 80

def get_axis_reps(axis_class, check = True,  verbose = 0):
    """Compute the global variables required in this module"""
    if check or verbose:
        VV = MMV(15)
        V_0 = Axis(1).in_space(VV)
    reps = axis_class.representatives()
    length = len(reps)
    indices = np.zeros((length,3), dtype = np.uint16)
    data = np.zeros((MAX_LEN_G + MAX_LEN_INVOLUTION) * length,
        dtype = np.uint32)
    prob =  np.zeros(length, dtype = np.double)
    start_involution = MM0(axis_class(1).g_axis)
    data_pos = 0
    d_sizes = axis_class.orbit_sizes()
    total_sizes = sum(d_sizes.values())
    size_factor = 1.0 / total_sizes
    sum_sizes = 0
    max_error = 1.0e-14 * total_sizes
    if verbose:
        print("Representatives of G_x0 orbits on axes")
    for i, (name, axis) in enumerate(reps.items()):
        g = chop_axis_representative(axis.g)
        assert len(g.mmdata) <= MAX_LEN_G 
        involution = aggressive_reduce_word(start_involution**MM0(g))
        assert len(involution.mmdata) <= MAX_LEN_INVOLUTION
        a = g.mmdata
        data[data_pos : data_pos + len(a)] = a
        indices[i,0] = data_pos 
        indices[i,1] = data_pos = data_pos + len(a)
        data_pos += len(a)
        a = involution.mmdata
        data[data_pos : data_pos + len(a)] = a
        indices[i,2] = data_pos = data_pos + len(a)
        data_pos += len(a)
        sum_sizes += d_sizes[name]
        prob[i] = sum_sizes * size_factor
        if check or verbose:
            g1 = MM('a', data[indices[i,0] : indices[i,1]])
            involution1 =  MM('a', data[indices[i,1] : indices[i,2]])
            assert MM(start_involution) ** MM(g1) == involution1
            assert V_0 * g1 ==  Axis(g).in_space(VV)
            p =  prob[i] - (prob[i - 1] if i else 0)
            assert abs(p * total_sizes - d_sizes[name]) < max_error
        if verbose:
            print("Orbit %s, relative size: %.3g" %
                 (name, PROB[i]))
            print("involution:\n", involution1)
    data = np.copy(data[:data_pos])
    return indices, data, prob



def get_vector_reps(axis_class, indices, data):
    V_0 = axis_class(1).in_space(VV)
    vectors = []
    for entry in indices:
        g = MM0('a', data[entry[0] : entry[1]])
        vectors.append(V_0 * g)
    return vectors




class RandomAxis:
    """Models a random 2A axis

    By default, the constructor of this class generates a uniform
    distributed random random 2A axis in the Monster. The functions
    in this class are optimized for utmost speed. Essentially, we
    compute the requested information about the axis on demand only.

    Each axis in this class has a number that describes the axis
    uniquely. If this number is passed to the constructor as a
    single argument, then the axis with that number is created.
    Method ``number`` returns the number of an axis. 
    """
    axis_class = Axis
    MAX_N_G_x0 = MAX_NUM_G_x0
    v15_start_name = V_AXIS_OPP 
    involution_start_name = G_AXIS_OPP
    indices, data, prob =  None, None, None
    _vectors = None
    n_orbits =  None

    __slots__ = 'orbit', 'a', 'num', '_v15'

    @classmethod
    def _get_vectors(cls):
        cls.indices, cls.data, cls.prob =  get_axis_reps(cls.axis_class)
        cls._vectors = get_vector_reps(cls.axis_class,
            cls.indices, cls.data)
        cls.n_orbits =  len(cls.indices)

    @staticmethod
    def _orbit_transformer(num = None):
        return rand_G_x0(num)

    @property
    def max_number(self):
        """Return maximum number of an axis minus 1.

        So for a axis ``a`` we have ``0 <= a.number < a.max_number``,
        for the constant value ``a.max_number``.
        """
        return self.MAX_N_G_x0 * self.n_orbits

    def __init__(self, number = None):
        if self.indices is None:
            self._get_vectors()
        if number is None:
            num = None
            self.orbit = gen_rng_uniform_to_dist(
                uniform(0, 1), self.n_orbits - 1, self.prob)
        else:   
            num, self.orbit = divmod(number, self.n_orbits) 
        self.a, self.num = self._orbit_transformer(num)
        self._v15 = None

    @property
    def number(self):
        """Return number of an axis.

        So for a axis ``a`` we have ``RandomAxis(a.number) == a``.
        So an axis genrated with this class can be indentified by its
        number. This number has no meaning outside class 
        ``RandomAxis``.

        Warning: The number returned by this method is implementation
        dependent!
        """
        return self.num * self.n_orbits + self.orbit

    def g(self, group = MM):
        r"""Group element mapping the start axis to the current axis 

        The function returns ``g`` with ``a.v15 = a.v15_start * a.g``.

        By default, ``g`` is an instance of class ``mmgroup.MM``. If
        parameter ``group`` is a different class implementing elements
        of the Monster (e.g class ``mmgroup.MM0``) then an instance of
        that class is returned. If ``group`` is the string ``'a'`` then
        the internal representation of an element of the Monster is
        returned as a numpy array of unsigned 32-bit integers.
        """
        ind = self.indices[self.orbit]
        orbit_rep = self.data[ind[0] : ind[1]]
        a = np.concatenate((orbit_rep, self.a))
        if group == 'a':
            return a
        return group('a', a)
        
    def g_axis(self, group = MM):
        r"""Return the 2A involution corresponding to the 2A axis

        The function returns the 2A involution ``g`` in the Monster
        corrsponding the axis.

        By default, ``g`` is an instance of class ``mmgroup.MM``. This
        may be changed by setting the optional parameter ``group``,
        which is interpreted as in method ``g``.
        """
        ind = self.indices[self.orbit]
        invol = self.data[ind[1] : ind[2]]
        a_inv = np.copy(self.a)
        mm_group_invert_word(a_inv, len(a_inv))
        a = np.concatenate((a_inv, invol, self.a))
        if group == 'a':
            return a
        return group('a', a)

    def g_axis_start(self, group = MM):
        r"""Return the fixed 2A involution :math:`\beta^+`.

        By default, the involution is an instance of class
        ``mmgroup.MM``. This may be changed by setting the optional
        parameter ``group``,  which is interpreted as in method ``g``.
        """
        if group == 'a':
            return MM(self.involution_start_name).mmdata
        return group(self.involution_start_name)


    def v15_start(cls):
        r"""Return the start axis vector of this class (mod 15)

        The method returns the fixed axis :math:`v^+`, where
        :math:`v^+` is the axis of the fixed involution
        :math:`\beta^+`.
        """
        return MMV(15)(cls.v15_start_name)

    def v15(self):
        r"""Return axis as vector with entries taken mod 15.

        For an instance ``a`` of this class the method returns the
        axis as a vector in the representation :math:`\rho` of the
        Monster with entries taken mod 15.
        """
        if self._v15:
            return self._v15
        if self._vectors is None:
            self._get_vectors()  
        self._v15 = v = self._vectors[self.orbit].copy()
        work =  mm_vector(v.p)
        mm_op_word(v.p, v.data, self.a, len(self.a), 1, work)
        return self._v15
        
    def axis(self):
        """Convert axis to instance of class ``mmgroup.axes.Axis``""" 
        return self.axis_class(self.g())



class RandomBabyAxis(RandomAxis):
    r"""Models a random 2A axis orthogonal to the standard axis

    By default, the constructor of this class generates a uniform
    distributed random random 2A axis in the Monster that is
    orthogonal to the standard axis :math:`v^+` corresponding to
    the standard 2A involution  :math:`\beta^+` in the Monster.

    The functions in this class are optimized for utmost speed.
    Essentially, we compute the requested information about the
    axis on demand only.

    Each axis in this class has a number that describes the axis
    uniquely. If this number is passed to the constructor as a
    single argument, then the axis with that number is created.
    Method ``number`` returns the number of an axis.
    The numbers obtained in that way are not compatible to
    the numbers computed in class ``RandomAxis``!

    The methods in this class are the same as in class
    ``RandomAxis``, except for the following differences:

    Method ``axis`` converts the current axis to instance of class
    ``mmgroup.axes.BabyAxis``

    Method ``g`` returns an element of the subgroup :math:`H^+` of
    the Monster mapping the start axis :math:`v^-` to the current
    axis. Here :math:`H^+` is the centralizer of the 2A involution
    :math:`\beta^+` (which is of structure :math:`2.B`).

    Method ``g_axis_start`` returns the fixed 2A involution
    :math:`\beta^-`.

    Method ``v15_start`` returns the start axis :math:`v^-` of
    this class (mod 15), which is the axis corresponding to the
    2A involution :math:`\beta^-`.
    """
    axis_class = BabyAxis
    MAX_N_G_x0 = MAX_NUM_H
    v15_start_name = V_AXIS_OPP 
    involution_start_name = G_AXIS_OPP
    indices, data, prob =  None, None, None
    _vectors = None
    n_orbits =  None

    @staticmethod
    def _orbit_transformer(num = None):
        return rand_G_x0_baby(num)

#############################################################
# Random group elements
#############################################################

def rand_mm_element(group_name = "M", group = MM):
    r"""Fast generation of a uniform random element of the Monster

    This function is not intended for public use!

    It is a low-level function returning a truly uniform
    distributed random element of the Monster.
    For some subgroups of the Monster such random elements can
    also be generated by specifying a string in parameter
    ``group_name`` as follows:

    =========== =======================================================================
    Value       Group
    =========== =======================================================================
    'M'         The Monster :math:`\mathbb{M}` (default)
    'B'         Centralizer :math:`2.B` of the standard 2A involution :math:`\beta^+`
    'G_x0'      Centralizer :math:`G_{x0}` of the standard 2B involution :math:`x_{-1}`
    'H'         The intersection :math:`2.B \cap G_{x0}`
    =========== =======================================================================


    By default, an instance of class ``mmgroup.MM`` is generated. If
    parameter ``group`` is a different class implementing elements
    of the Monster (e.g class ``mmgroup.MM0``) then an instance of
    that class is returned. If ``group`` is the string ``'a'`` then
    the internal representation of an element of the Monster is
    returned as a numpy array of unsigned 32-bit integers.
    """
    if group_name == "G_x0":
        mm =  rand_G_x0()[0]
    elif group_name in ["H", "B", "M"]:
        mm = rand_G_x0_baby()[0]
        if group_name in "BM":
            mm1 = RandomBabyAxis().g('a')
            mm = np.concatenate((mm, mm1))
        if group_name == "M":
            mm2 =  RandomAxis().g('a')
            mm = np.concatenate((mm, mm2))
    else:
        ERR = "Illegal name %s for subgroup of Monster"
        raise ValueError(ERR % group_name)
    if group == 'a':
        return mm
    return group('a', mm)
        


