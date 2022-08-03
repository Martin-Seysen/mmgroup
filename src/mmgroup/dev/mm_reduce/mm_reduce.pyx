
# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


from numbers import Integral
import numpy as np
from libc.stdint cimport uint64_t, uint32_t, int32_t, uint8_t, int8_t
from libc.stdlib cimport malloc, free

from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
from mmgroup.structures.construct_mm import iter_strings_from_atoms



cimport cython
include "../pxd_files/mm_reduce.pxi"
cimport generators as gen
cimport mm_reduce as mr

ctypedef gt_word_type *p_gt_word_type
ctypedef gt_subword_type *p_gt_subword_type


cdef uint32_t[:] mm_as_array_view(a):
    cdef uint32_t[:] a_view 
    cdef uint32_t length 
    if isinstance(a, AbstractMMGroupWord):
        a1 = a.mmdata
        assert isinstance(a1, np.ndarray)
    else:
        a1 = np.array(a, dtype = np.uint32)
    a_view = a1
    return a_view



cdef class GtWord():
    cdef p_gt_word_type p_gt 
    group_name = "GtWord"

    @cython.boundscheck(False)    
    def __cinit__(self, arg):
        cdef uint32_t[:] a_view
        cdef uint32_t length
        if isinstance(arg, Integral):
            length = max(arg, 1)
        else:
            a_view = mm_as_array_view(arg)
            length = len(a_view)
            length = mr.gt_word_n_subwords(&a_view[0], length)
        self.p_gt = <p_gt_word_type>malloc(mr.gt_word_size(length)) 
        mr.gt_word_init(self.p_gt, length)


    def  __dealloc__(self):
        if self.p_gt != NULL: (self.p_gt)

    @cython.boundscheck(False)    
    def __init__(self, arg):
        cdef int32_t res
        cdef uint32_t[:] a_view
        cdef uint32_t length
        if self.p_gt == NULL:
            err = "Memory allocation has failed in class GtWord"
            raise MemoryError(err)
        if not isinstance(arg, Integral):
            a_view = mm_as_array_view(arg)
            length = len(a_view)
            res = mr.gt_word_load(self.p_gt, &a_view[0], length)
            assert res == 0

    def fst(self):
        cdef int i = self.p_fst - self.a
        return GtSubWord(self, i)

    def end(self):
        cdef uint32_t i = self.p_end - self.a
        return GtSubWord(self, i)

    def subwords(self):
        cdef uint32_t i
        cdef p_gt_subword_type cur = self.p_gt.p_fst
        lst = [self.fst]
        while cur.p_next != NULL and cur.type == 1:
            cur = cur.p_next
            i = cur - self.p_gt.a
            lst.append(GtSubWord(self, i))
        return lst

    def size(self):
        return self.p_gt.n_subwords

    @cython.boundscheck(False)    
    def mmdata(self, group = None):
        cdef int32_t length = mr.gt_word_length(self.p_gt)
        assert length >= 0
        a = np.zeros(length, dtype = np.uint32)
        cdef uint32_t[:] a_view = a
        mr.gt_word_store(self.p_gt, &a_view[0], length)
        if group is None:
           return a
        assert issubclass(group, AbstractMMGroupWord)
        return group('a', a)

    def __str__(self):
        """Represent group element as a string"""
        strings = iter_strings_from_atoms(self.mmdata, abort_if_error=0)
        s = "*".join(strings) if len(strings) else "1"
        return "%s<%s>" % (self.group_name, s)




cdef class GtSubWord(object):
    cdef p_gt_subword_type p_gtsub
    cdef GtWord master
    group_name = "GtSubWord"
 
    def __init__(self, GtWord master, int i):
        self.master = master
        assert 0 <= i < master.n_subwords
        cdef p_gt_subword_type p0 = self.master.p_gt.a 
        self.p_gtsub = p0 + i

    def next(self):
        cdef uint32_t t = self.p_gtsub.type 
        assert t == 1
        cdef p_gt_subword_type p_next = self.p_gtsub.p_next
        cdef p_gt_subword_type p0 = self.master.p_gt.a 
        cdef int i = p_next - p0
        return GtSubWord(self.master, i) 

    def prev(self):
        cdef uint32_t t = self.p_gtsub.type 
        assert t > 0 
        cdef p_gt_subword_type p_prev = self.p_gtsub.p_prev
        cdef p_gt_subword_type p0 = self.master.p_gt.a 
        cdef int i = p_prev - p0
        return GtSubWord(self.master, i) 

    def data(self):
        cdef uint32_t *data = &(self.p_gtsub.data)[0]  
        cdef uint32_t length = self.p_gtsub.length 
        a = np.empty(length, dtype = np.uint32)
        cdef int j
        for j in range(length): a[j]  = data[j]
        return a

    def isendmark(self):
        cdef uint32_t t = self.p_gtsub.type 
        assert 1 <= t <= 2
        return bool(t - 1)

    def img_Omega(self, inverse = False):
        cdef int inv = bool(inverse)
        cdef uint32_t length = self.p_gtsub.length
        cdef int32_t res = gen.gen_leech2_op_word_leech2(0x800000, 
             self.p_gtsub.data, length, inv)
        assert 0 < res < 0x1000000, res
        if inv:
            self.p_gtsub.img_Omega_inv = res
        else:
            self.p_gtsub.img_Omega = res
        return res


    def setdata(self, data):
        cdef int j, length = len(data)
        assert length < mr.MAX_GT_WORD_DATA
        for j in range(length): self.p_gtsub.data[j]  = data[j]

    def delete(self):
        assert self.type == 1
        cdef p_gt_subword_type p = mr.gt_subword_delete(self.p_gtsub) 
        cdef p_gt_subword_type p0 = self.master.p_gt.a 
        cdef int i = p - p0
        return GtSubWord(self.master, i) 
     
    def insert(self):
        assert self.type > 0
        cdef p_gt_subword_type p = mr.gt_subword_insert(self.p_gtsub) 
        if (p == NULL):
            err = "GtSubWord object is out of local memeory"
            raise ValueError(err)
        cdef p_gt_subword_type p0 = self.master.p_gt.a 
        cdef int i = p - p0
        return GtSubWord(self.master, i)

    @cython.boundscheck(False)    
    def load(self, arg):
        cdef uint32_t[:] a_view = mm_as_array_view(arg)
        cdef uint32_t length = len(a_view)
        cdef int32_t res = mr.gt_subword_load(
             self.p_gtsub, &a_view[0], length)
        assert res >= 0
        return self

    @cython.boundscheck(False)    
    def mmdata(self, group = None):
        if group is None:
           return self.data
        assert issubclass(group, AbstractMMGroupWord)
        return group('a', self.data)

    def __str__(self):
        """Represent group element as a string"""
        strings = iter_strings_from_atoms(self.mmdata, abort_if_error=0)
        s = "*".join(strings) if len(strings) else "1"
        return "%s<%s>" % (self.group_name, s)


    def __eq__(self, GtSubWord other):
        cdef p_gt_subword_type p_me = self.p_gtsub
        cdef p_gt_subword_type p_other = other.p_gtsub
        return p_me == p_other and p_me.p_master == p_other.p_master

 
    def __ne__(self, GtSubWord other):
        return self.__eq__(other)



