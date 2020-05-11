# This .pxd file has been generated automatically. Do not edit!


# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from libc.stdint cimport uint32_t, uint16_t, uint8_t, int32_t
from libc.stdint cimport uint16_t
from libc.stdint cimport uint64_t as uint_mmv_t



cdef extern from "mm_op3.h":
    void mm_op3_pi(uint_mmv_t *v_in, uint32_t delta, uint32_t pi, uint_mmv_t *v_out)
    void mm_op3_delta(uint_mmv_t *v_in, uint32_t delta, uint_mmv_t *v_out)
    uint32_t mm_op3_copy(uint_mmv_t *mv1, uint_mmv_t *mv2)
    uint32_t mm_op3_compare(uint_mmv_t *mv1, uint_mmv_t *mv2)
    void mm_op3_vector_add(uint_mmv_t *mv1, uint_mmv_t *mv2)
    void mm_op3_scalar_mul(int32_t factor, uint_mmv_t *mv1)
    void mm_op3_xy(uint_mmv_t *v_in, uint32_t f, uint32_t e, uint32_t eps, uint_mmv_t *v_out)
    void mm_op3_t(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
    void mm_op3_xi(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
    void mm_op3_group_n(uint_mmv_t *v, uint32_t *g, uint_mmv_t *work)
    void mm_op3_word(uint_mmv_t *v, uint32_t *g, int32_t len_g, int32_t e, uint_mmv_t *work)
