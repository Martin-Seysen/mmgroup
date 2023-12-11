# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from libc.stdint cimport uint64_t, uint32_t, uint16_t, uint8_t



cdef extern from "mat24_functions.h":
    const uint32_t MAT24_BASIS[24]


import numpy


"""yet to be documented!!!! """


cimport cython
cimport generators as g



include "generators.pxi"


def make_table(uint32_t u_box, uint32_t u_exp):
    assert 1 <= u_box <= 5 and 1 <= u_exp <= 2
    a = numpy.zeros(32768, dtype = numpy.uint16) 
    cdef uint16_t[::1] pa = a
    cdef uint32_t length
    length = g.gen_xi_make_table(u_box, u_exp, &pa[0])
    return a[:length]


def invert_table(table, uint32_t n_columns, uint32_t len_result):
    result = numpy.zeros(len_result, dtype = numpy.uint16) 
    cdef uint32_t length = len(table)
    assert length & 31 == 0 and len_result & 31 == 0
    cdef uint16_t[::1] pr = result
    cdef uint16_t[::1] pt = table
    g.gen_xi_invert_table(&pt[0], length, n_columns, &pr[0], len_result)
    return result


def split_table(table, modulus):
    cdef uint32_t length = len(table)
    assert length & 31 == 0
    a = numpy.zeros(length >> 5, dtype = numpy.uint32) 
    cdef uint32_t[::1] pa = a
    cdef uint16_t[::1] pt = table
    g.gen_xi_split_table(&pt[0], length, modulus, &pa[0])
    return table, a


#####################################################################################
# The random generator
#####################################################################################

from numbers import Integral
import numpy as np

ctypedef uint64_t* p_uint64_t


g.gen_rng_seed_init()

try:
    # The following statement requires python version >= 3.8
    from threading import get_native_id
except ImportError:
    # This works for python version >= 3.3
    from threading import get_ident
    from os import getpid
    def get_native_id():
        return (getpid() << 32) +  get_ident()





_seed_dict = {}


cdef _make_std_rand_gen():
    global _seed_dict
    cdef uint64_t thread_id = get_native_id()
    a = np.zeros(4, dtype = np.uint64)
    cdef uint64_t[::1] a1 = a
    g.gen_rng_seed(&a1[0])
    _seed_dict[thread_id] = a
    return a



@cython.wraparound(False)
@cython.boundscheck(False)
cdef p_uint64_t c_rng_get_seed(seed):
    cdef uint64_t[::1] a1
    if seed is None:
        try:
            a1 = _seed_dict[get_native_id()]
            return &a1[0]
        except KeyError:
            a1 = _make_std_rand_gen()
            return &a1[0]
    else:
        a1 = seed
        assert len(seed) >= 4
        return &a1[0]
       
def rand_get_seed():
    cdef uint64_t[::1] a1
    try:
        return _seed_dict[get_native_id()]
    except KeyError:
        return _make_std_rand_gen()




def rand_make_seed(value):
    a = np.zeros(4, dtype = np.uint64)
    cdef uint64_t[::1] a1 = a
    g.gen_rng_seed_no(&a1[0], value)
    return a
        


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef rand_bytes_modp(uint32_t p, uint32_t num_bytes, seed = None):
    cdef uint32_t p1 = p
    cdef uint32_t n1 = num_bytes
    a = np.zeros(n1, dtype = np.uint8)
    cdef uint8_t[::1] a1 = a 
    cdef p_uint64_t p_seed = c_rng_get_seed(seed)
    cdef int32_t status
    status = g.gen_rng_bytes_modp(p1, &a1[0], n1, p_seed)
    if (status):
        err = "Bad modulus for generating random bytes"
        raise ValueError(err)
    return a

    
@cython.wraparound(False)
@cython.boundscheck(False)
cpdef rand_fill_bytes_modp(p, array_bytes, seed = None):
    cdef uint32_t p1 = p
    cdef uint32_t n1 = len(array_bytes)
    cdef uint8_t[::1] a1 = array_bytes 
    cdef p_uint64_t p_seed = c_rng_get_seed(seed)
    cdef int32_t status
    status = g.gen_rng_bytes_modp(p1, &a1[0], n1, p_seed)
    if (status):
        err = "Bad modulus for generating random bytes"
        raise ValueError(err)
    
cpdef uint32_t rand_gen_modp(uint32_t p, seed = None):
    assert 1 <= p <= 0x100000000
    cdef uint32_t p1 = p & 0xffffffff
    cdef p_uint64_t p_seed = c_rng_get_seed(seed)
    cdef uint32_t result = g.gen_rng_modp(p1, p_seed)
    return result


cpdef uint64_t rand_gen_bitfields_modp(p, uint32_t d, seed = None):
    assert 1 <= p <= 0x10000000000000000
    cdef uint32_t p1 = p & 0xffffffffffffffff
    cdef p_uint64_t p_seed = c_rng_get_seed(seed)
    cdef uint64_t result = g.gen_rng_bitfields_modp(p1, d, p_seed)
    return result




     