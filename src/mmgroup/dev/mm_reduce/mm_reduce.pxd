from libc.stdint cimport uint64_t, uint32_t, uint16_t, uint8_t
from libc.stdint cimport int64_t, int32_t
from libc.stdint cimport uint%{INT_BITS}_t as uint_mmv_t

INT_BITS = %{INT_BITS}


cdef extern from "mm_reduce.h":
    enum: MAX_GT_WORD_DATA
    enum: MM_COMPRESS_TYPE_NENTRIES

    ctypedef struct gt_subword_type:
        uint32_t eof  
        uint32_t length
        uint32_t img_Omega
        uint32_t t_exp
        uint32_t reduced
        gt_subword_type *p_prev
        gt_subword_type *p_next
        uint32_t data[MAX_GT_WORD_DATA]

    ctypedef struct gt_word_type:
        gt_subword_type *p_end
        gt_subword_type *p_node
        gt_subword_type *p_free
        int32_t reduce_mode;

    ctypedef struct mm_compress_type:
        uint64_t nx
        uint32_t w[MM_COMPRESS_TYPE_NENTRIES]




