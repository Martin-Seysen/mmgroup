# This .pxd file has been generated automatically. Do not edit!



from libc.stdint cimport uint32_t, uint16_t, uint8_t, int32_t
from libc.stdint cimport uint16_t
from libc.stdint cimport uint64_t as uint_mmv_t

INT_BITS = 64

cdef extern from "mm_basics.h":
    void mm_aux_read_mmv1(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len)
    void mm_aux_read_direct_mmv1(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len)
    void mm_aux_write_mmv1(uint32_t p, uint8_t *b, uint_mmv_t *mv, uint32_t len)
    void mm_aux_read_mmv24(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len)
    void mm_aux_write_mmv24(uint32_t p, uint8_t *b, uint_mmv_t *mv, uint32_t len)
    uint8_t mm_aux_get_mmv1(uint32_t p, uint_mmv_t *mv, uint32_t i)
    void mm_aux_put_mmv1(uint32_t p, uint8_t value, uint_mmv_t *mv, uint32_t i)
    uint32_t mm_aux_mmv_size(uint32_t p)
    void mm_aux_zero_mmv(uint32_t p, uint_mmv_t *mv)
    uint8_t mm_aux_get_mmv(uint32_t p, uint_mmv_t *mv, uint32_t i)
    void mm_aux_put_mmv(uint32_t p, uint8_t value, uint_mmv_t *mv, uint32_t i)
    void mm_aux_random_mmv(uint32_t p, uint8_t *seed, uint_mmv_t *mv)
    int32_t mm_aux_reduce_mmv(uint32_t p, uint_mmv_t *mv)
    int32_t mm_aux_check_mmv(uint32_t p, uint_mmv_t *mv)
    void mm_aux_small24_expand(uint8_t *b_src, uint8_t *b_dest)
    void mm_aux_small24_compress(uint8_t *b_src, uint8_t *b_dest)
    void mm_aux_mmv_to_bytes(uint32_t p, uint_mmv_t *mv, uint8_t *b)
    void mm_aux_bytes_to_mmv(uint32_t p, uint8_t *b, uint_mmv_t *mv)
    int32_t mm_aux_mmv_to_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp)
    void mm_aux_mmv_extract_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp, uint32_t length)
    void mm_aux_mmv_add_sparse(uint32_t p, uint32_t *sp, uint32_t length, uint_mmv_t *mv)
    void mm_aux_mmv_set_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp, uint32_t length)
    uint32_t mm_aux_index_extern_to_sparse(uint32_t i)
    void mm_aux_array_extern_to_sparse(uint32_t *a, uint32_t len)
    int32_t mm_aux_index_sparse_to_extern(uint32_t i)
    uint32_t mm_sparse_purge(uint32_t *sp, uint32_t length)
    uint32_t mm_sparse_bitsort(uint32_t *a, uint32_t len, uint32_t mask)
    uint32_t mm_sparse_alloc_reduce(uint32_t p, uint32_t *sp, uint32_t length)
    uint32_t mm_sparse_reduce(uint32_t p, uint32_t *sp, uint32_t length)
    void mm_sparse_mul(uint32_t p, uint32_t *sp, uint32_t length, int32_t factor)
    void mm_rng_seed(uint8_t *seed, uint32_t no, uint8_t *key, uint32_t len)
    void mm_rng_gen_modp(uint8_t p, uint8_t *seed, uint8_t *out, uint32_t len)
    void mm_rng_gen_bytes(uint8_t *seed, uint8_t *out, uint32_t len)
    void  mm_sub_test_prep_pi_64(uint32_t delta, uint32_t pi, uint32_t *p_tbl)
    void  mm_sub_test_prep_xy(uint32_t f, uint32_t e, uint32_t eps, uint32_t n, uint32_t *p_tbl)
    void mm_group_n_mul_delta_pi(uint32_t *g, uint32_t delta, uint32_t pi)
    void mm_group_n_mul_inv_delta_pi(uint32_t *g, uint32_t delta, uint32_t pi)
    void mm_group_n_mul_x(uint32_t * g, uint32_t e)
    void mm_group_n_mul_y(uint32_t * g, uint32_t f)
    void mm_group_n_mul_t(uint32_t * g, uint32_t exp)
    void mm_group_n_mul_element(uint32_t *g, uint32_t *g1)
    void mm_group_n_mul_inv_element(uint32_t *g, uint32_t *g1)
    void mm_group_n_clear(uint32_t *g)
    uint32_t mm_group_n_mul_atom(uint32_t *g, uint32_t atom)
    uint32_t mm_group_n_to_word(uint32_t *g, uint32_t *word)
    uint32_t mm_group_split_word_n(uint32_t *word, uint32_t length, uint32_t *g)
    uint32_t mm_group_mul_words(uint32_t *w1, uint32_t l1, uint32_t *w2, uint32_t l2, int32_t e)
    uint32_t mm_sub_get_table_xi(uint32_t e, uint32_t i, uint32_t j, uint32_t k)
