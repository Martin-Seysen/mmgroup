

from libc.stdint cimport uint32_t, uint16_t, uint8_t, int32_t
from libc.stdint cimport uint16_t
from libc.stdint cimport uint64_t as uint_mmv_t

INT_BITS = 64


######################################################################
### Wrappers for C functions from file mm_basics.pxd
######################################################################


cimport mm_basics


def mm_aux_read_mmv1(p, mv, b, len):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint8_t[::1] b_v_ = b
    cdef uint32_t len_v_ = len
    mm_basics.mm_aux_read_mmv1(p_v_, &mv_v_[0], &b_v_[0], len_v_)

def mm_aux_read_direct_mmv1(p, mv, b, len):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint8_t[::1] b_v_ = b
    cdef uint32_t len_v_ = len
    mm_basics.mm_aux_read_direct_mmv1(p_v_, &mv_v_[0], &b_v_[0], len_v_)

def mm_aux_write_mmv1(p, b, mv, len):
    cdef uint32_t p_v_ = p
    cdef uint8_t[::1] b_v_ = b
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t len_v_ = len
    mm_basics.mm_aux_write_mmv1(p_v_, &b_v_[0], &mv_v_[0], len_v_)

def mm_aux_read_mmv24(p, mv, b, len):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint8_t[::1] b_v_ = b
    cdef uint32_t len_v_ = len
    mm_basics.mm_aux_read_mmv24(p_v_, &mv_v_[0], &b_v_[0], len_v_)

def mm_aux_write_mmv24(p, b, mv, len):
    cdef uint32_t p_v_ = p
    cdef uint8_t[::1] b_v_ = b
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t len_v_ = len
    mm_basics.mm_aux_write_mmv24(p_v_, &b_v_[0], &mv_v_[0], len_v_)

def mm_aux_get_mmv1(p, mv, i):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t i_v_ = i
    cdef uint8_t ret_
    ret_ = mm_basics.mm_aux_get_mmv1(p_v_, &mv_v_[0], i_v_)
    return ret_

def mm_aux_put_mmv1(p, value, mv, i):
    cdef uint32_t p_v_ = p
    cdef uint8_t value_v_ = value
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t i_v_ = i
    mm_basics.mm_aux_put_mmv1(p_v_, value_v_, &mv_v_[0], i_v_)

def mm_aux_mmv_size(p):
    cdef uint32_t p_v_ = p
    cdef uint32_t ret_
    ret_ = mm_basics.mm_aux_mmv_size(p_v_)
    return ret_

def mm_aux_zero_mmv(p, mv):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    mm_basics.mm_aux_zero_mmv(p_v_, &mv_v_[0])

def mm_aux_get_mmv(p, mv, i):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t i_v_ = i
    cdef uint8_t ret_
    ret_ = mm_basics.mm_aux_get_mmv(p_v_, &mv_v_[0], i_v_)
    return ret_

def mm_aux_put_mmv(p, value, mv, i):
    cdef uint32_t p_v_ = p
    cdef uint8_t value_v_ = value
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t i_v_ = i
    mm_basics.mm_aux_put_mmv(p_v_, value_v_, &mv_v_[0], i_v_)

def mm_aux_random_mmv(p, seed, mv):
    cdef uint32_t p_v_ = p
    cdef uint8_t[::1] seed_v_ = seed
    cdef uint_mmv_t[::1] mv_v_ = mv
    mm_basics.mm_aux_random_mmv(p_v_, &seed_v_[0], &mv_v_[0])

def mm_aux_reduce_mmv(p, mv):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef int32_t ret_
    ret_ = mm_basics.mm_aux_reduce_mmv(p_v_, &mv_v_[0])
    return ret_

def mm_aux_check_mmv(p, mv):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef int32_t ret_
    ret_ = mm_basics.mm_aux_check_mmv(p_v_, &mv_v_[0])
    return ret_

def mm_aux_small24_expand(b_src, b_dest):
    cdef uint8_t[::1] b_src_v_ = b_src
    cdef uint8_t[::1] b_dest_v_ = b_dest
    mm_basics.mm_aux_small24_expand(&b_src_v_[0], &b_dest_v_[0])

def mm_aux_small24_compress(b_src, b_dest):
    cdef uint8_t[::1] b_src_v_ = b_src
    cdef uint8_t[::1] b_dest_v_ = b_dest
    mm_basics.mm_aux_small24_compress(&b_src_v_[0], &b_dest_v_[0])

def mm_aux_mmv_to_bytes(p, mv, b):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint8_t[::1] b_v_ = b
    mm_basics.mm_aux_mmv_to_bytes(p_v_, &mv_v_[0], &b_v_[0])

def mm_aux_bytes_to_mmv(p, b, mv):
    cdef uint32_t p_v_ = p
    cdef uint8_t[::1] b_v_ = b
    cdef uint_mmv_t[::1] mv_v_ = mv
    mm_basics.mm_aux_bytes_to_mmv(p_v_, &b_v_[0], &mv_v_[0])

def mm_aux_mmv_to_sparse(p, mv, sp):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t[::1] sp_v_ = sp
    cdef int32_t ret_
    ret_ = mm_basics.mm_aux_mmv_to_sparse(p_v_, &mv_v_[0], &sp_v_[0])
    return ret_

def mm_aux_mmv_extract_sparse(p, mv, sp, length):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t[::1] sp_v_ = sp
    cdef uint32_t length_v_ = length
    mm_basics.mm_aux_mmv_extract_sparse(p_v_, &mv_v_[0], &sp_v_[0], length_v_)

def mm_aux_mmv_add_sparse(p, sp, length, mv):
    cdef uint32_t p_v_ = p
    cdef uint32_t[::1] sp_v_ = sp
    cdef uint32_t length_v_ = length
    cdef uint_mmv_t[::1] mv_v_ = mv
    mm_basics.mm_aux_mmv_add_sparse(p_v_, &sp_v_[0], length_v_, &mv_v_[0])

def mm_aux_mmv_set_sparse(p, mv, sp, length):
    cdef uint32_t p_v_ = p
    cdef uint_mmv_t[::1] mv_v_ = mv
    cdef uint32_t[::1] sp_v_ = sp
    cdef uint32_t length_v_ = length
    mm_basics.mm_aux_mmv_set_sparse(p_v_, &mv_v_[0], &sp_v_[0], length_v_)

def mm_aux_index_extern_to_sparse(i):
    cdef uint32_t i_v_ = i
    cdef uint32_t ret_
    ret_ = mm_basics.mm_aux_index_extern_to_sparse(i_v_)
    return ret_

def mm_aux_array_extern_to_sparse(a, len):
    cdef uint32_t[::1] a_v_ = a
    cdef uint32_t len_v_ = len
    mm_basics.mm_aux_array_extern_to_sparse(&a_v_[0], len_v_)

def mm_aux_index_sparse_to_extern(i):
    cdef uint32_t i_v_ = i
    cdef int32_t ret_
    ret_ = mm_basics.mm_aux_index_sparse_to_extern(i_v_)
    return ret_

def mm_sparse_purge(sp, length):
    cdef uint32_t[::1] sp_v_ = sp
    cdef uint32_t length_v_ = length
    cdef uint32_t ret_
    ret_ = mm_basics.mm_sparse_purge(&sp_v_[0], length_v_)
    return ret_

def mm_sparse_bitsort(a, len, mask):
    cdef uint32_t[::1] a_v_ = a
    cdef uint32_t len_v_ = len
    cdef uint32_t mask_v_ = mask
    cdef uint32_t ret_
    ret_ = mm_basics.mm_sparse_bitsort(&a_v_[0], len_v_, mask_v_)
    return ret_

def mm_sparse_alloc_reduce(p, sp, length):
    cdef uint32_t p_v_ = p
    cdef uint32_t[::1] sp_v_ = sp
    cdef uint32_t length_v_ = length
    cdef uint32_t ret_
    ret_ = mm_basics.mm_sparse_alloc_reduce(p_v_, &sp_v_[0], length_v_)
    return ret_

def mm_sparse_reduce(p, sp, length):
    cdef uint32_t p_v_ = p
    cdef uint32_t[::1] sp_v_ = sp
    cdef uint32_t length_v_ = length
    cdef uint32_t ret_
    ret_ = mm_basics.mm_sparse_reduce(p_v_, &sp_v_[0], length_v_)
    return ret_

def mm_sparse_mul(p, sp, length, factor):
    cdef uint32_t p_v_ = p
    cdef uint32_t[::1] sp_v_ = sp
    cdef uint32_t length_v_ = length
    cdef int32_t factor_v_ = factor
    mm_basics.mm_sparse_mul(p_v_, &sp_v_[0], length_v_, factor_v_)

def mm_rng_seed(seed, no, key, len):
    cdef uint8_t[::1] seed_v_ = seed
    cdef uint32_t no_v_ = no
    cdef uint8_t[::1] key_v_ = key
    cdef uint32_t len_v_ = len
    mm_basics.mm_rng_seed(&seed_v_[0], no_v_, &key_v_[0], len_v_)

def mm_rng_gen_modp(p, seed, out, len):
    cdef uint8_t p_v_ = p
    cdef uint8_t[::1] seed_v_ = seed
    cdef uint8_t[::1] out_v_ = out
    cdef uint32_t len_v_ = len
    mm_basics.mm_rng_gen_modp(p_v_, &seed_v_[0], &out_v_[0], len_v_)

def mm_rng_gen_bytes(seed, out, len):
    cdef uint8_t[::1] seed_v_ = seed
    cdef uint8_t[::1] out_v_ = out
    cdef uint32_t len_v_ = len
    mm_basics.mm_rng_gen_bytes(&seed_v_[0], &out_v_[0], len_v_)

def mm_sub_test_prep_pi_64(delta, pi, p_tbl):
    cdef uint32_t delta_v_ = delta
    cdef uint32_t pi_v_ = pi
    cdef uint32_t[::1] p_tbl_v_ = p_tbl
    mm_basics.mm_sub_test_prep_pi_64(delta_v_, pi_v_, &p_tbl_v_[0])

def mm_sub_test_prep_xy(f, e, eps, n, p_tbl):
    cdef uint32_t f_v_ = f
    cdef uint32_t e_v_ = e
    cdef uint32_t eps_v_ = eps
    cdef uint32_t n_v_ = n
    cdef uint32_t[::1] p_tbl_v_ = p_tbl
    mm_basics.mm_sub_test_prep_xy(f_v_, e_v_, eps_v_, n_v_, &p_tbl_v_[0])

def mm_group_n_mul_delta_pi(g, delta, pi):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t delta_v_ = delta
    cdef uint32_t pi_v_ = pi
    mm_basics.mm_group_n_mul_delta_pi(&g_v_[0], delta_v_, pi_v_)

def mm_group_n_mul_inv_delta_pi(g, delta, pi):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t delta_v_ = delta
    cdef uint32_t pi_v_ = pi
    mm_basics.mm_group_n_mul_inv_delta_pi(&g_v_[0], delta_v_, pi_v_)

def mm_group_n_mul_x(g, e):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t e_v_ = e
    mm_basics.mm_group_n_mul_x(&g_v_[0], e_v_)

def mm_group_n_mul_y(g, f):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t f_v_ = f
    mm_basics.mm_group_n_mul_y(&g_v_[0], f_v_)

def mm_group_n_mul_t(g, exp):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t exp_v_ = exp
    mm_basics.mm_group_n_mul_t(&g_v_[0], exp_v_)

def mm_group_n_mul_element(g, g1):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t[::1] g1_v_ = g1
    mm_basics.mm_group_n_mul_element(&g_v_[0], &g1_v_[0])

def mm_group_n_mul_inv_element(g, g1):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t[::1] g1_v_ = g1
    mm_basics.mm_group_n_mul_inv_element(&g_v_[0], &g1_v_[0])

def mm_group_n_clear(g):
    cdef uint32_t[::1] g_v_ = g
    mm_basics.mm_group_n_clear(&g_v_[0])

def mm_group_n_mul_atom(g, atom):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t atom_v_ = atom
    cdef uint32_t ret_
    ret_ = mm_basics.mm_group_n_mul_atom(&g_v_[0], atom_v_)
    return ret_

def mm_group_n_to_word(g, word):
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t[::1] word_v_ = word
    cdef uint32_t ret_
    ret_ = mm_basics.mm_group_n_to_word(&g_v_[0], &word_v_[0])
    return ret_

def mm_group_split_word_n(word, length, g):
    cdef uint32_t[::1] word_v_ = word
    cdef uint32_t length_v_ = length
    cdef uint32_t[::1] g_v_ = g
    cdef uint32_t ret_
    ret_ = mm_basics.mm_group_split_word_n(&word_v_[0], length_v_, &g_v_[0])
    return ret_

def mm_group_mul_words(w1, l1, w2, l2, e):
    cdef uint32_t[::1] w1_v_ = w1
    cdef uint32_t l1_v_ = l1
    cdef uint32_t[::1] w2_v_ = w2
    cdef uint32_t l2_v_ = l2
    cdef int32_t e_v_ = e
    cdef uint32_t ret_
    ret_ = mm_basics.mm_group_mul_words(&w1_v_[0], l1_v_, &w2_v_[0], l2_v_, e_v_)
    return ret_

def mm_sub_get_table_xi(e, i, j, k):
    cdef uint32_t e_v_ = e
    cdef uint32_t i_v_ = i
    cdef uint32_t j_v_ = j
    cdef uint32_t k_v_ = k
    cdef uint32_t ret_
    ret_ = mm_basics.mm_sub_get_table_xi(e_v_, i_v_, j_v_, k_v_)
    return ret_

