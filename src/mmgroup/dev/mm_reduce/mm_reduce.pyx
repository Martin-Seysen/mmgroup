
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
include "mm_reduce.pxi"
cimport generators as gen
cimport mm_reduce as mr

ctypedef gt_word_type *p_gt_word_type
# ctypedef gt_subword_type *p_gt_subword_type


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

    @staticmethod
    def _complain(res, method):
        if (res < 0):
            err = "Internal error %s in class GtWord, method %s"
            raise ValueError(err % (hex(res), method))  

    @cython.boundscheck(False)    
    def __cinit__(self, *args, **kwds):
        self.p_gt = mr.gt_word_alloc(0, NULL, 0)
        if self.p_gt == NULL:
            err = "Memory allocation has failed in class GtWord"
            raise MemoryError(err)


    def  __dealloc__(self):
        mr.gt_word_free(self.p_gt)

    @cython.boundscheck(False)    
    def __init__(self, arg = None, reduce_mode = 1):
        cdef int32_t res
        cdef uint32_t[:] a_view
        cdef uint32_t length
        self.p_gt.reduce_mode = reduce_mode
        if arg is not None:
            a_view = mm_as_array_view(arg)
            length = len(a_view)
            res = mr.gt_word_append(self.p_gt, &a_view[0], length)
            if res < 0:
                self._complain(res, "__init__")

    def eof(self):
        return self.p_gt.p_node.eof

    def reduce_mode(self):
        return self.p_gt.reduce_mode

    def set_reduce_mode(self, uint32_t reduce_mode):
        self.p_gt.reduce_mode = reduce_mode

    def seek(self, int32_t pos, int32_t seek_set):
        cdef int32_t res = mr.gt_word_seek(self.p_gt, pos, seek_set)
        if res < 0:
            self._complain(res, "seek")



    @cython.boundscheck(False)    
    def mmdata(self, group = None):
        cdef int32_t length = mr.gt_word_length(self.p_gt)
        assert length >= 0
        a = np.zeros(length, dtype = np.uint32)
        cdef uint32_t[:] a_view = a
        length = mr.gt_word_store(self.p_gt, &a_view[0], length)
        if (length < 0):
            self._complain(length, "mmdata")
        a = a[:length]      
        if group is None:
            return a
        assert issubclass(group, AbstractMMGroupWord)
        return group('a', a)


    def rule_join(self):
        cdef int32_t res = mr.gt_word_rule_join(self.p_gt)
        if res < 0:
            self._complain(res, "rule_join")
        return res

    def rule_t_xi_t(self):
        cdef int32_t res = mr.gt_word_rule_t_xi_t(self.p_gt)
        if res < 0:
            s = "Error %d in method GtWord.rule_t_xi_t:" % res
            self.display_subwords(s)
            self._complain(res, "rule_join")
        return res

    @cython.boundscheck(False)    
    def append_sub_part(self, a):
        cdef uint32_t[:] a_view = mm_as_array_view(a)
        cdef int32_t n = mr.gt_word_append_sub_part(
            self.p_gt, &a_view[0], len(a_view))
        if n < 0:
            self._complain(n, "append_sub_part")
        return n, tail_word(a, n)


    @cython.boundscheck(False)    
    def append(self, a):
        cdef uint32_t[:] a_view = mm_as_array_view(a)
        cdef int32_t n = mr.gt_word_append(
            self.p_gt, &a_view[0], len(a_view))
        if n < 0:
            self._complain(n, "append")

    def reduce_sub(self, uint32_t mode = 1):
        cdef int32_t res = mr.gt_word_reduce_sub(self.p_gt, mode)
        if res < 0:
            self._complain(res, "reduce_sub")

    def reduce(self):
        cdef int32_t res = mr.gt_word_reduce(self.p_gt)
        if res < 0:
            self._complain(res, "reduce_sub")

    @cython.boundscheck(False)    
    def subwords(self):
        cdef int32_t maxlength = mr.gt_word_len_pickle(self.p_gt)
        assert maxlength > 0
        buf = np.zeros(maxlength, dtype = np.uint32)
        cdef uint32_t[:] b_view = buf
        cdef int32_t length = mr.gt_word_pickle(
            self.p_gt, &b_view[0], maxlength)
        assert length > 0
        cdef uint32_t n = b_view[1], i, a_len, tmp
        cdef int32_t fpos = -1
        subword_list = []
        for i in range(n):
            a = buf[buf[i+2]:]
            sub = GtSubWord()
            if a[0]:
                fpos = i
            sub.img_Omega = tmp = a[1]
            sub.t_exp = tmp = a[2]
            sub.reduced = tmp = a[3]
            sub.length = a_len = a[4]
            sub.data = a[5 : 5 + a_len]
            subword_list.append(sub)
        return fpos, subword_list
  
    @cython.boundscheck(False)
    def as_int(self):
        self.reduce()
        cdef uint64_t a[4]
        cdef uint64_t[:] a_view = a
        cdef res = mr.gt_word_compress(self.p_gt, &a_view[0])
        if res:
             raise ValueError("Cound not convert element to int")
        return (int(a[0]) + (int(a[1]) << 64) +
             (int(a[2]) << 128) + (int(a[3]) << 192))

    @cython.boundscheck(False)
    def as_int_debug_compress(self):
        cdef mm_compress_type c
        a = np.zeros(MM_COMPRESS_TYPE_NENTRIES + 1, dtype = np.uint64)
        cdef int32_t status = mr.gt_word_to_mm_compress(self.p_gt, &c)
        cdef uint32_t i, j=1, ww
        if (status < 0):
             err = "Function gt_word_to_mm_compress failed, status = %d"
             raise ValueError(err % status)
        a[0] = c.nx
        for i in range(MM_COMPRESS_TYPE_NENTRIES):
            ww = c.w[i]
            if ww:
                a[j] =  ww
                j += 1
        return a[:j]


    def __len__(self):
        cdef int32_t l = mr.gt_word_length(self.p_gt)
        if l < 0:
            self._complain(l, "__len__")
        return l


    def __str__(self):
        """Represent group element as a string"""
        strings = iter_strings_from_atoms(self.mmdata(), abort_if_error=0)
        s = "*".join(strings)
        if len(s) == 0: strings = "1"
        return "%s<%s>" % (self.group_name, s)


    @cython.boundscheck(False)    
    def _insert(self):
        cdef int32_t res = mr.gt_word_insert(self.p_gt)
        if res < 0:
            s = "Error %d in method GtWord._insert:" % res
            self.display_subwords(s)
            self._complain(res, "_insert")

    @cython.boundscheck(False)    
    def _delete(self):
        cdef int32_t res = mr.gt_word_delete(self.p_gt)
        if res < 0:
            s = "Error %d in method GtWord._delete:" % res
            self.display_subwords(s)
            self._complain(res, "delete")


    def display_subwords(self, text = None):
        if text:
            print(text)
        fpos, subword_list = self.subwords()
        for i, w in enumerate(subword_list):
            m = "^" if i == fpos else " "
            m += "r " if w.reduced else "  "
            strings = iter_strings_from_atoms(w.data, abort_if_error=0)
            s = "*".join(strings)
            if w.t_exp:
                s += ("*" if len(s) else "") + "t^%d" % w.t_exp
            print(m + s)
        m = "^end" if fpos == -1 else " end"
        print(m)


        



