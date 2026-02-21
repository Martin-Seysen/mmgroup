from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import warnings
import copy
import numbers
import re
import numpy as np
from numbers import Integral
from random import randint


#from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup import MMV, MMVector, Xsp2_Co1
from mmgroup.mm_op import INT_BITS, PROTECT_OVERFLOW

from mmgroup.generators import gen_leech3_neg, gen_leech3_reduce
from mmgroup.clifford12 import error_string, chk_qstate12
from mmgroup.clifford12 import xsp2co1_rep_mod3_conv_mm_op
from mmgroup.clifford12 import xsp2co1_rep_mod3_unit_vector
from mmgroup.clifford12 import xsp2co1_rep_mod3_mul_word
from mmgroup.clifford12 import xsp2co1_rep_mod3_mul_elem
from mmgroup.clifford12 import xsp2co1_to_vect_mod3
from mmgroup.structures.xsp2_co1 import get_error_pool
from mmgroup.structures.xsp2_co1 import str_xsp2_co1
from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.abstract_rep_space import AbstractRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractRepSpace
from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord

MMSpace3 =  MMV(3)



######################################################################
# Creating a vector in rep 4096x as an instance of class QStateMatrix
######################################################################

ERR_QSTATE_SHAPE = "QStateMatrix object must have shape (0,12) here"

def obj_to_leech_mod3(obj, unit):
    if obj == 'r':
        obj = [randint(0, 2) for i in range(24)]
    elif isinstance(obj, dict):
        ll = [0] * 24
        for key, value in obj.items():
            ll[key] = value
        obj = ll
    if isinstance(obj, Integral):
        if unit:
            assert 0 <= obj < 24
            return 1 << obj
        else:
            return  obj & 0xffffffffffff
    elif isinstance(obj, list):
        v3 = 0
        for i, x in enumerate(obj[:24]):
            x = x % 3
            if x:
                v3 |= 1 << (i + 24 * (x - 1))
        return v3
    else:
        ERR = "Cannot construct Leech lattice vector from a %s object."
        raise ValueError(ERR % type(tag))

def _err_missing_args(self):
    ERR = "Missing argument in constructor of class Xsp2_Co1_Vector"
    raise ValueError(ERR)


def create_vector(data, *args):
    # data must be a numpy array of shape = (1,), dtype = np.uint64,
    # and length 14.
    np.copyto(data, 0)
    if len(args) == 0 or not args[0]:
         return data
    factor = 1
    if isinstance(args[0], Integral):
        factor = args[0] % 3
        if factor == 0:
            return data
        args = args[1:]
    if len(args) == 0:
         _err_missing_args()
    tag = args[0]
    if isinstance(tag, str) and tag in "ZY":
        if len(args) < 3:
            _err_missing_args()
        if isinstance(args[1], Integral):
            i0 = args[1] + ((tag == "Y") << 11)
        elif args[1] == 'r':
            i0 = randint(0, 0xfff)
        else:
            ERR = "Unit vector index for Xsp2_Co1_Vector must be integral"
            raise TypeError(ERR)
        v3 = obj_to_leech_mod3(args[2], unit=True)
        xsp2co1_rep_mod3_unit_vector(i0 & 0xfff, v3, data)
    elif isinstance(tag, QStateMatrix):
        if len(args) < 2:
            _err_missing_args()
        v3 = obj_to_leech_mod3(args[1], unit=False)
        a = tag._xsp2co1_rep_mod3_from_qs(v3)
        np.copyto(data, a)
    else:
        ERR = "Cannot construct an Xsp2_Co1_Vector from a %s object."
        raise TypeError(ERR % type(tag))
    if factor == 2:
        data[0] = gen_leech3_neg(data[0])
    return data


######################################################################
# Modelling a vector
######################################################################

class Xsp2_Co1_Vector(AbstractMmRepVector):
    r"""Model a vector in the space :math:`4096_x \otimes 24_x` (mod 3)

    This lass can represent only tensors :math:`q \otimes v_3` (mod 3)
    where :math:`q` is a quadratic state vector in :math:`4096_x`,
    and :math:`v_3` is a vector in the Leech lattice mod 3. Note that
    :math:`4096_x \otimes 24_x` (mod 3) is a subspace of the
    representation and :math:`\rho_3` of the monster. Vectors in
    this class are stored in a much more compact form than vectors
    in representation :math:`\rho_3`.

    The constructor of this class can be used to create unit vectors
    in :math:`4096_x \otimes 24_x` (mod 3) using the same language as
    for the construction of vectors in :math:`\rho_3` in class ``MMV``.
    So, calling

    ``Xsp2_Co1_Vector(tag, i, j)`` ,

    where tag is 'X' or 'Y' contructs the same unit vector as
    calling  ``MMV(3, tag, i, j)``, for 0 <= i < 2048, 0 <= j < 24.
    But in this class there is a wraparound so that
    ``Xsp2_Co1_Vector('Y', i, j)`` is equivalent to
    ``Xsp2_Co1_Vector('X', i + 2048, j)``. Any tags different
    from 'X' and 'Y' represent unit vectors in :math:`\rho_3` that
    are not in :math:`4096_x \otimes 24_x`; so they are illegal here.
 
    We may also construct a tensor product in this class from an
    element ``q`` of :math:`4096_x` and a vector ``v3`` in
    :math:`24_x` by calling

    ``Xsp2_Co1_Vector(q, v3)`` .

    Then ``q`` must be an instance of class ``QStateMatrix``
    representing a quadratic state vector of shape (0, 12) or
    (12, 0). Parameter ``v3`` must be a vector in the Leech
    lattice mod 3 in *Leech lattice mod 3 encoing'*

    The first argument ``tag`` or ``q`` in the constructor may be
    preceded by an (optional) integer argument, which will be
    interpreted as a scalar factor, if present.

    More options for parameter ``v3`` are yet t0 be documented!
    """
    __slots__ = '_data', 'space', '_mmv'
    def __init__(self, space):
        __slots__ = 'space', '_data'
        self.space = space
        self._data = np.zeros(14, dtype = np.uint64)
        self._mmv = None

    def _set_vector(self, *args):
        create_vector(self._data, *args)
        self._mmv = None

    @property
    def mmv(self):
        if self._mmv is None:
            self._mmv = MMSpace3()
            START = 116416//32
            chk_qstate12(xsp2co1_rep_mod3_conv_mm_op(self._data,
                 self._mmv.data[START:]))
        return self._mmv

    @property
    def is_zero(self):
        v3 = int(self._data[0])
        return not bool((v3 ^ (v3 >> 24)) & 0xffffff)

    def check(self):
        """Check if the vector is correct

        Raise ValueError if the vector is errorneous.
        """
        self.space.check(self)

    def dump(self):
        self.space.dump_vector(self)    

    def as_sparse(self):
        """Yet to be documented!!

        """
        vm = self.mmv
        return vm.space().as_sparse(vm)

    def as_tuples(self, v1):
        r"""Return vector in tuple representation

        The function returns a list of tuples ``(factor, tag, i0, i1)``.
        Here ``(tag, i0, i1)`` describes a basis vector and ``value``
        is the coordinate of that vector.

        Entries with coordinate zero are dropped.
        """
        return self.mmv.as_tuples()

    def mul_Gx0(self, g):
        r"""Alternative implementation of multiplcation with group element

        The function replaces the the tensor ``v`` in this object
        by ``v * g`` for an element ``g`` of :math:`G_{x0}`. This
        function works only if the factor of the tensor product ``v``
        in the space  :math:`24_x` is a *short* vector in the Leech
        lattice (mod 3), i.e. a vector of type 2.
        """
        g = Xsp2_Co1(g)
        chk_qstate12(xsp2co1_rep_mod3_mul_elem(self._data, g._data))
        return v1

    def factors(self):
        if self.is_zero:
            return  QStateMatrix(), 0
        aq = self._data[1:14]
        n = len(aq)
        while n > 0 and aq[n-1] == 0:
            n = n - 1
        qs = QStateMatrix(0, 12, aq[:n], 0)
        v3 = gen_leech3_reduce(int(self._data[0]))
        return qs, v3

    def show(self):
        if self.is_zero:
            return "<Xsp2_Co1_Vector 0>\n"
        qs, v3 = self.factors()
        v3r = xsp2co1_to_vect_mod3(v3)
        s = "<Xsp2_Co1_Vector [%s] (x)\n" % str_mmv3(v3r)
        s += qs.show(1, 0) + ">"
        return s

    def __str__(self):
        return self.show()


######################################################################
# class SparseMmSpace
######################################################################




class Xsp2_Co1_Space(AbstractRepSpace):
    """Models the vector space for class Xsp2_Co1_Vector

    YET TO BE DOCUMENTED !!!

    """

    vector_type = Xsp2_Co1_Vector
    mmspace = MMSpace3

    def __init__(self):
        """Create a 196884-dimensional representation of the monster

        All calculations are done modulo the odd number p
        """
        super(Xsp2_Co1_Space, self).__init__(3, Xsp2_Co1)


    ### Standard methods that must be overwritten ####################

    def zero(self):
        return self.vector_type(self)

    def unit(self, *data):
        v = self.vector_type(self)
        v._set_vector(*data)
        return v

    def iadd(self, v1, v2):
        if v1.is_zero:
            return v2
        if v2.is_zero:
            return v1
        ERR = "Vector addition not supported in this space"
        raise TypeError(ERR)

    def imul_scalar(self, v1, a):
        a = int(a % 3)
        if a == 0:
            v1._mmv = None
            np.copyto(self._data, 0)
        elif a == 2:
            v1._mmv = None
            v1._data[0] = gen_leech3_neg(v1._data[0])
        return v1

    def imul_group_word(self, v1, g):
        if v1.is_zero:
            return v1
        v1._mmv = None
        assert isinstance(g, AbstractMMGroupWord)
        a = g.mmdata
        chk_qstate12(xsp2co1_rep_mod3_mul_word(v1._data, a, len(a)))
        return v1

    def vector_get_item(self, v1, index):
        return v1.mmv[index]

    def vector_set_item(self, v1, index, value) :
        ERR = "Item assigment not supported in this vector space"
        raise TypeError(ERR)

    def equal_vectors(self, v1, v2):
        if v1.is_zero != v2.is_zero:
            return False
        if v1.is_zero:
            return True
        self.reduce(v1)
        self.reduce(v2)
        return (v1._data == v2._data).all()


    ### Standard methods that must be overwritten ####################

    def copy_vector(self, v1):
        """Return deep copy of group element v1"""
        assert v1.space == self
        v2 = self.zero()
        v2._data[:] = v1._data
        return v2

    def reduce(self, v):
        if v.is_zero:
            np.copyto(v._data, 0)
        return v

    def str_vector(self, v):
        """Convert vector v to a string

        """
        if v.is_zero:
            return "Xsp2 vector 0"
        return "Xsp2 vector " + str_xsp2_co1(v._data, factor=64, t=True)

    ### Standard methods that need not be overwritten #################

    def parse(self, s):
        ERR = "Cannot convert a string to a  Xsp2_Co1_Vector"
        raise TypeError(ERR)

    ### Extra methods ##################################################

    def dump_vector(self, v1):
        a = np.zeros(4096, dtype = np.uint64)    
        chk_qstate12(xsp2co1_rep_mod3_conv_mm_op(self._data, a))
        or_sum = 0
        print("Xsp2_Co1 vector:")
        for i, x in enumerate(a):
            x = int(x) & 0xffffffffffff
            if (x): 
                print("0x%04x: 0x%012x" % (i, x))
            or_sum |= x
        if not or_sum:
            print("Vector is equal to 0")

    def check(self, v1):
        pass  # Yet nothing do do here!!
        return v1

Space_ZY =  Xsp2_Co1_Space()


######################################################################
# Display a vector in SparseMmSpace
######################################################################

def str_mmv3(v3):
    v3 = int(v3)
    lst = []
    for i in range(3):
        s = "".join([".+-."[(v3 >> j) & 3] for j in range(0, 16, 2)])
        lst.append(s)
        v3 >>= 16
    return " ".join(lst)

def dump_zy(text, v):
    print(text + " (hex):")
    if isinstance(v, Xsp2_Co1_Vector):
        v = v.mmv
    assert isinstance(v, MMVector), type(v)
    assert v.p == 3
    v.reduce()
    START = 116416//32
    zy = v.data[START : START + 4096]
    n = 0
    for i, v3 in enumerate(zy):
        if v3:
            e = "\n" if n & 1 else "  "
            print(" %3x: %s" % (i, str_mmv3(v3)), end = e)
            n += 1
    if (n == 0):
        print(" (zero vector)", end = "")
        n += 1
    if (n & 1):
        print()

