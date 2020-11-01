# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from libc.stdint cimport uint32_t, uint16_t, uint8_t

#from cpython cimport array as c_array



include "../pxd_files/mat24_xi_functions.pxd"


cdef extern from "mat24_functions.h":
    const uint32_t MAT24_BASIS[24]


import numpy


"""yet to be documented!!!! """



def xi_g_gray(uint32_t v1):
    return mat24_xi_g_gray(v1)


def xi_w2_gray(uint32_t v1):
    return mat24_xi_w2_gray(v1)


def xi_g_cocode(uint32_t v1):
    return mat24_xi_g_cocode(v1)


def xi_w2_cocode(uint32_t v1):
    return mat24_xi_w2_cocode(v1)


def xi_op_xi(uint32_t x1, uint32_t u_exp):
    return mat24_xi_op_xi(x1, u_exp)


def xi_mul_leech(uint32_t x1, uint32_t x2):
    return mat24_xi_mul_leech(x1, x2)


def xi_leech_to_short(uint32_t x1):
    return mat24_xi_leech_to_short(x1)


def xi_short_to_leech(uint32_t x1):
    return mat24_xi_short_to_leech(x1)


def xi_op_xi_short(uint32_t x1, uint32_t u_exp):
    return mat24_xi_op_xi_short(x1, u_exp)


def make_table(uint32_t u_box, uint32_t u_exp):
    assert 1 <= u_box <= 5 and 1 <= u_exp <= 2
    a = numpy.zeros(32768, dtype = numpy.uint16) 
    cdef uint16_t[::1] pa = a
    cdef uint32_t length
    length = mat24_xi_make_table(u_box, u_exp, &pa[0])
    return a[:length]


def invert_table(table, uint32_t n_columns, uint32_t len_result):
    result = numpy.zeros(len_result, dtype = numpy.uint16) 
    cdef uint32_t length = len(table)
    assert length & 31 == 0 and len_result & 31 == 0
    cdef uint16_t[::1] pr = result
    cdef uint16_t[::1] pt = table
    mat24_xi_invert_table(&pt[0], length, n_columns, &pr[0], len_result)
    return result


def split_table(table, modulus):
    cdef uint32_t length = len(table)
    assert length & 31 == 0
    a = numpy.zeros(length >> 5, dtype = numpy.uint32) 
    cdef uint32_t[::1] pa = a
    cdef uint16_t[::1] pt = table
    mat24_xi_split_table(&pt[0], length, modulus, &pa[0])
    return table, a



