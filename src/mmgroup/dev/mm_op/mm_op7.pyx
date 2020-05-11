
# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from libc.stdint cimport uint32_t, uint16_t, uint8_t, int32_t
from libc.stdint cimport uint16_t
from libc.stdint cimport uint64_t as uint_mmv_t




######################################################################
### Wrappers for C functions from file mm_op7.pxd
######################################################################


cimport mm_op7


def op_pi(v_in, delta, pi, v_out):
    cdef uint_mmv_t[::1] v_in_v_ = v_in
    cdef uint32_t delta_v_ = delta
    cdef uint32_t pi_v_ = pi
    cdef uint_mmv_t[::1] v_out_v_ = v_out
    mm_op7.mm_op7_pi(&v_in_v_[0], delta_v_, pi_v_, &v_out_v_[0])

def op_delta(v_in, delta, v_out):
    cdef uint_mmv_t[::1] v_in_v_ = v_in
    cdef uint32_t delta_v_ = delta
    cdef uint_mmv_t[::1] v_out_v_ = v_out
    mm_op7.mm_op7_delta(&v_in_v_[0], delta_v_, &v_out_v_[0])

def op_copy(mv1, mv2):
    cdef uint_mmv_t[::1] mv1_v_ = mv1
    cdef uint_mmv_t[::1] mv2_v_ = mv2
    cdef uint32_t ret_
    ret_ = mm_op7.mm_op7_copy(&mv1_v_[0], &mv2_v_[0])
    return ret_

def op_compare(mv1, mv2):
    cdef uint_mmv_t[::1] mv1_v_ = mv1
    cdef uint_mmv_t[::1] mv2_v_ = mv2
    cdef uint32_t ret_
    ret_ = mm_op7.mm_op7_compare(&mv1_v_[0], &mv2_v_[0])
    return ret_

def op_vector_add(mv1, mv2):
    cdef uint_mmv_t[::1] mv1_v_ = mv1
    cdef uint_mmv_t[::1] mv2_v_ = mv2
    mm_op7.mm_op7_vector_add(&mv1_v_[0], &mv2_v_[0])

def op_scalar_mul(factor, mv1):
    cdef int32_t factor_v_ = factor
    cdef uint_mmv_t[::1] mv1_v_ = mv1
    mm_op7.mm_op7_scalar_mul(factor_v_, &mv1_v_[0])

def op_xy(v_in, f, e, eps, v_out):
    cdef uint_mmv_t[::1] v_in_v_ = v_in
    cdef uint32_t f_v_ = f
    cdef uint32_t e_v_ = e
    cdef uint32_t eps_v_ = eps
    cdef uint_mmv_t[::1] v_out_v_ = v_out
    mm_op7.mm_op7_xy(&v_in_v_[0], f_v_, e_v_, eps_v_, &v_out_v_[0])

def op_t(v_in, exp, v_out):
    cdef uint_mmv_t[::1] v_in_v_ = v_in
    cdef uint32_t exp_v_ = exp
    cdef uint_mmv_t[::1] v_out_v_ = v_out
    mm_op7.mm_op7_t(&v_in_v_[0], exp_v_, &v_out_v_[0])

def op_xi(v_in, exp, v_out):
    cdef uint_mmv_t[::1] v_in_v_ = v_in
    cdef uint32_t exp_v_ = exp
    cdef uint_mmv_t[::1] v_out_v_ = v_out
    mm_op7.mm_op7_xi(&v_in_v_[0], exp_v_, &v_out_v_[0])

def op_group_n(v, g, work):
    cdef uint_mmv_t[::1] v_v_ = v
    cdef uint32_t[::1] g_v_ = g
    cdef uint_mmv_t[::1] work_v_ = work
    mm_op7.mm_op7_group_n(&v_v_[0], &g_v_[0], &work_v_[0])

def op_word(v, g, len_g, e, work):
    cdef uint_mmv_t[::1] v_v_ = v
    cdef uint32_t[::1] g_v_ = g
    cdef int32_t len_g_v_ = len_g
    cdef int32_t e_v_ = e
    cdef uint_mmv_t[::1] work_v_ = work
    mm_op7.mm_op7_word(&v_v_[0], &g_v_[0], len_g_v_, e_v_, &work_v_[0])


######################################################################
### Constants
######################################################################


MMV_ENTRIES = 247488

INT_BITS = 64

LOG_INT_BITS = 6

P = 7

FIELD_BITS = 4

LOG_FIELD_BITS = 2

INT_FIELDS = 16

LOG_INT_FIELDS = 4

P_BITS = 3

MMV_INTS = 15468

