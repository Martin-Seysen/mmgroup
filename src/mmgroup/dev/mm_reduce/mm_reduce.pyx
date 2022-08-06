
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



def tail_word(a, n):
    if isinstance(a, AbstractMMGroupWord):
        a1 = a.mmdata[n:]
        assert isinstance(a1, np.ndarray)
        return a.__class__('a', a1)
    else:
        return np.array(a, dtype = np.uint32)[n:]

@cython.boundscheck(False)    
def _n_subwords(w):
    cdef uint32_t[:] a_view = mm_as_array_view(w)
    cdef uint32_t  n = len(a_view)
    n = mr.gt_word_n_subwords(&a_view[0], n)
    return n


class GtSubWord():
    """Container for data in subword of instance of class GtWord"""
    pass



cdef class GtWord():
    cdef p_gt_word_type p_gt 
    group_name = "GtWord"

    @staticmethod
    def n_subwords(w):
        return _n_subwords(w)

    @cython.boundscheck(False)    
    def __cinit__(self, arg):
        cdef uint32_t length
        if isinstance(arg, Integral):
            length = max(arg, 1)
        else:
            length = _n_subwords(arg) + 6
        self.p_gt = <p_gt_word_type>malloc(mr.gt_word_size(length)) 
        if self.p_gt == NULL:
            err = "Memory allocation has failed in class GtWord"
            raise MemoryError(err)
        mr.gt_word_init(self.p_gt, length)


    def  __dealloc__(self):
        if self.p_gt != NULL: free(self.p_gt)

    @cython.boundscheck(False)    
    def __init__(self, arg):
        cdef int32_t res
        cdef uint32_t[:] a_view
        cdef uint32_t length
        if not isinstance(arg, Integral):
            a_view = mm_as_array_view(arg)
            length = len(a_view)
            res = mr.gt_word_append(self.p_gt, &a_view[0], length)
            assert res == 0

    def eof(self):
        return self.p_gt.p_node.eof

    def seek(self, int32_t pos, int32_t seek_set):
        cdef int32_t res = mr.gt_word_seek(self.p_gt, pos, seek_set)
        if res < 0:
            err = "Method seek failed on <GtWord>  object"
            raise ValueError(err) 


    def size(self):
        return self.p_gt.n_subwords

    @cython.boundscheck(False)    
    def mmdata(self, group = None):
        cdef int32_t length = mr.gt_word_length(self.p_gt)
        assert length >= 0
        a = np.zeros(length, dtype = np.uint32)
        cdef uint32_t[:] a_view = a
        length = mr.gt_word_store(self.p_gt, &a_view[0], length)
        if (length < 0):
            err = "Method mmdata failed on <GtWord>  object"
            raise ValueError(err)
        a = a[:length]      
        if group is None:
            return a
        assert issubclass(group, AbstractMMGroupWord)
        return group('a', a)

    @cython.boundscheck(False)    
    def mmdata_sub(self, group = None):
        cdef int32_t length = mr.MAX_GT_WORD_DATA
        a = np.zeros(length, dtype = np.uint32)
        cdef uint32_t[:] a_view = a
        length = mr.gt_word_store_sub(self.p_gt, &a_view[0], length)
        if (length < 0):
            err = "Method mmdata failed on <GtWord>  object"
            raise ValueError(err)
        a = a[:length]      
        if group is None:
           return a
        assert issubclass(group, AbstractMMGroupWord)
        return group('a', a)

    def join(self):
        cdef int32_t res = mr.gt_word_join(self.p_gt)
        assert res >= 0

    @cython.boundscheck(False)    
    def append_sub_part(self, a):
        cdef uint32_t[:] a_view = mm_as_array_view(a)
        cdef int32_t n = mr.gt_word_append_sub_part(
            self.p_gt, &a_view[0], len(a_view))
        assert n >= 0
        return n, tail_word(a, n)

    @cython.boundscheck(False)    
    def append_sub(self, a):
        cdef uint32_t[:] a_view = mm_as_array_view(a)
        cdef int32_t n = mr.gt_word_append_sub(
            self.p_gt, &a_view[0], len(a_view))
        assert n >= 0

    @cython.boundscheck(False)    
    def append(self, a):
        cdef uint32_t[:] a_view = mm_as_array_view(a)
        cdef int32_t n = mr.gt_word_append(
            self.p_gt, &a_view[0], len(a_view))
        assert n >= 0

    def reduce_sub(self, uint32_t mode = 1):
        cdef int32_t res = mr.gt_word_reduce_sub(self.p_gt, mode)
        assert res >= 0

    def reduce(self, uint32_t mode = 1):
        cdef int32_t res = mr.gt_word_reduce(self.p_gt, mode)
        assert res >= 0


    @cython.boundscheck(False)    
    def subwords(self):
        cdef p_gt_subword_type ptr = self.p_gt.p_end.p_next
        cdef uint32_t i
        cdef uint32_t[:] a_view
        cdef uint32_t[:] src_view
        cdef uint32_t length
        cdef int32_t maxlength = mr.MAX_GT_WORD_DATA
        subword_list = []
        while not ptr.eof:
            sub = GtSubWord()
            sub.data = np.zeros(maxlength, dtype = np.uint32)
            a_view = sub.data
            src_view = ptr.data
            length = ptr.length
            for i in range(length):
                a_view[i] = src_view[i]
            sub.data = sub.data[:length]
            sub.length = length
            sub.t_exp = ptr.t_exp
            sub.img_Omega = ptr.img_Omega
            sub.reduced = ptr.reduced
            subword_list.append(sub)
            ptr = ptr.p_next
        return subword_list

    
    def __len__(self):
        cdef int32_t l = mr.gt_word_length(self.p_gt)
        return l


    def __str__(self):
        """Represent group element as a string"""
        strings = iter_strings_from_atoms(self.mmdata(), abort_if_error=0)
        s = "*".join(strings)
        if len(s) == 0: strings = "1"
        return "%s<%s>" % (self.group_name, s)






