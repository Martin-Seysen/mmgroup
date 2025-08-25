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


#from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup import MM0
from mmgroup.mm_space import MMSpace, MMV
from mmgroup.mm_op import INT_BITS, PROTECT_OVERFLOW

from mmgroup.generators import gen_leech2to3_short
from mmgroup.clifford12 import xsp2co1_elem_to_qs, xsp2co1_qs_to_elem_i 
from mmgroup.clifford12 import error_string, chk_qstate12
from mmgroup.clifford12 import xsp2co1_neg_elem, xsp2co1_elem_row_mod3
from mmgroup.clifford12 import xsp2co1_mul_elem, xsp2co1_reduce_elem
from mmgroup.structures.xsp2_co1 import get_error_pool
from mmgroup.structures.xsp2_co1 import Xsp2_Co1, str_xsp2_co1
from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.abstract_rep_space import AbstractRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_rep_space import AbstractRepSpace
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
#from mmgroup.structures.mm_space_indices import sparse_to_tuples




MMSpace3 =  MMV(3)
MMGroup3 = MM0

######################################################################
# Creating a leech lattice vector mod 3 from a tuple
######################################################################

ERR_LEECH = "Tuple does not describe a short Leech lattice vector %s"


def tuple_to_leech_mod3(tag, i0 = -1, i1 = -1):
    leech_mod2 = MMSpace.index_to_short_mod2(tag, i0, i1)
    res = gen_leech2to3_short(leech_mod2)
    if (res == 0):
        err = get_error_pool(15)
        if len(err):
            print("tuple_to_leech_mod3 input:", hex(leech_mod2))
            print("Leech", MMSpace.index_to_short(tag, i0, i1))
        raise ValueError(ERR_LEECH % "mod 3")
    return res


CODE_MOD3 = {0:0, 1:1, 2:0x1000001}
def array_to_leech_mod3(a):
    return sum(CODE_MOD3[x % 3] << i for i, x in enumerate(a[:24]))


######################################################################
# Creating a vector in rep 4096x as an instance of class QStateMatrix
######################################################################

ERR_QSTATE_SHAPE = "QStateMatrix object must have shape (0,12) here"

def obj_to_qstate(obj):
    if isinstance(obj, QStateMatrix):
        if obj.shape == (0,12):
            return obj.copy().extend_zero(12,0).reshape((12,12))
        if obj.shape == (12,12):
            return obj
        raise ValueError(ERR_QSTATE_SHAPE)
    if isinstance(obj, Integral):
        res = QStateMatrix(0, 12,  obj & 0xfff)
        if (obj & 0x1000):
            res = - res 
        return res.extend_zero(0,12).reshape((12,12))
    ERR = "Illegal type for creating an Xsp2_Co1_Vector"
    raise TypeError(ERR)  


######################################################################
# Modelling a vector
######################################################################

class Xsp2_Co1_Vector(AbstractMmRepVector):
    def __init__(self, space):
        __slots__ = 'space', '_data', 'is_zero'
        self.space = space
        self.is_zero = True
        self._data = np.zeros(26, dtype = np.uint64)

    @property
    def short3(self):
        return int(self._data[0])

    @property
    def short3_list(self):
        x = self._data[0]
        return  [(((x >> i) & 1) + ((x >> (23 + i)) & 2)) % 3 
            for i in range(24)]

    @property
    def qs(self):
        return QStateMatrix(xsp2co1_elem_to_qs(self._data))

    def _set_vector(self, rep_4096, rep_leech):
        qs = obj_to_qstate(rep_4096)
        if isinstance(rep_leech, Integral): 
            x = rep_leech
        else:
            tag = rep_leech[0]
            if isinstance(tag, str):
                x = tuple_to_leech_mod3(*rep_leech)
            else:
                x = array_to_leech_mod3(rep_leech)
        self._data[:] = xsp2co1_qs_to_elem_i(qs, x)
        self.is_zero = False
       
    def _set_zero(self):
        self.is_zero = True

    def check(self):
        """Check if the vector is correct

        Raise ValueError if the vector is errorneous.
        """
        self.space.check(self)
        
    def as_tuples(self):
        return self.space.as_tuples(self)   
        
    def dump(self):
        self.space.dump_vector(self)    

    def as_mmspace_vector(self):
        return self.space.as_mmspace_vector(self) 



######################################################################
# class SparseMmSpace
######################################################################




class Xsp2_Co1_Space(AbstractRepSpace):
    """Models the sparse representation 198884x of the monster group. 

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

    def unit(self, rep_4096, rep_leech):
        v = self.vector_type(self)
        v._set_vector(rep_4096, rep_leech)
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
        v1.is_zero =  v1.is_zero or a == 0
        if a == 2 and not v1.is_zero:
            xsp2co1_neg_elem(v1._data)
        return v1

    def imul_group_word(self, v1, g):
        if v1.is_zero:
            return v1
        chk_qstate12(xsp2co1_mul_elem(v1._data, g._data, v1._data))
        return v1

    def vector_get_item(self, v1, index):
        return self.as_mmspace_vector(v1)[index]

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
        v2.is_zero = v1.is_zero
        return v2

    def reduce(self, v):
        if not v.is_zero:
            chk_qstate12(xsp2co1_reduce_elem(v._data))
        return v


    def str_vector(self, v):
        """Convert vector v to a string

        """
        if v.is_zero:
            return "Xsp2 vector 0"
        return "Xsp2 vector " + str_xsp2_co1(v._data, factor=64, t=True)

    ### Standard methods that need not be overwritten #################

    def parse(self, s):
        ERR = "Cannot convert a strig to a  Xsp2_Co1_Vector"
        raise TypeError(ERR)

    ### Extra methods ##################################################

    def as_mmspace_vector(self, v1):
        vm = MMV(3)()
        assert v1.qs.shape == (12,12) 
        if INT_BITS == 64:
            START = 116416//32
            assert len(vm.data) >=  START + 4096
            chk_qstate12(xsp2co1_elem_row_mod3(v1._data, 0, 
                 vm.data[START:]))
            return vm
        ERR = "Function not supported on %d-bit systems"
        raise NotImplementedError(ERR % INT_BITS)

    def as_sparse(self, v1):
        """Yet to be documented!!

        """
        vm = self.as_mmspace_vector(v1)
        return vm.space().as_sparse(vm)

    def as_tuples(self, v1):
        """Return vector in tuple representation

        The function returns a list of tuples ``(factor, tag, i0, i1)``.
        Here ``(tag, i0, i1)`` describes a basis vector and ``value``
        is the coordinate of that vector.

        Entries with coordinate zero are dropped.
        """
        vm = self.as_mmspace_vector(v1)
        return vm.as_tuples()
        
    def dump_vector(self, v1):
        a = np.zeros(4096, dtype = np.uint64)    
        chk_qstate12(xsp2co1_elem_row_mod3(v1._data, 0, a))
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
